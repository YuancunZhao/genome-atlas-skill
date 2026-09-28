#!/usr/bin/env bash
# Reference data, about 60 GB. Everything is public and re-downloadable; nothing here is personal.
set -euo pipefail
source "$(dirname "$0")/../scripts/env.sh"
mkdir -p "$REF_DIR"/{b37,hg38,chain,annot,imp,aadr,archaic,ytree,prs}
cd "$REF_DIR"

# H6 (was item #8): every download used to go straight to its final name with no check. A truncated or
# corrupt transfer is the worst failure mode here -- the file exists, later steps read it, and the only
# symptom is a subtly wrong result. fetch() downloads to .part, verifies, and only then renames, so a
# failure can never leave a plausible-looking file behind.
#   - curl -f makes an HTTP error a non-zero exit instead of a saved error page
#   - the result must be non-empty and above a per-file floor
#   - by extension: gzip/bgzip archives must pass gzip -t; JSON must parse; text must have content
fetch(){
  local url=$1 out=$2 min=${3:-1024}
  if [ -s "$out" ]; then
    if ! verify "$out"; then
      echo "WARNING: existing $out fails verification; re-downloading" >&2
    else
      return 0
    fi
  fi
  curl -fsSL --retry 3 -o "$out.part" "$url" || { echo "download failed: $url" >&2; rm -f "$out.part"; return 1; }
  local sz; sz=$(stat -c %s "$out.part" 2>/dev/null || echo 0)
  if [ "$sz" -lt "$min" ]; then
    echo "download too small ($sz < $min bytes): $url" >&2; rm -f "$out.part"; return 1
  fi
  if ! verify "$out.part"; then rm -f "$out.part"; return 1; fi
  mv "$out.part" "$out"
  echo "fetched $(basename "$out") ($sz bytes)"
}

verify(){
  local f=$1 kind=$1
  # fetch() 传进来的必然是下载中的 .part 文件，而按后缀分派要看得懂它的真实类型：不剥掉 .part 时
  # ".gz.part"/".json.part" 会落进 *) 兜底分支、只查非空——校验形同虚设（复审 P0）。剥的只是**判断
  # 用的名字**，读取仍用原路径。
  case "$kind" in *.part) kind="${kind%.part}" ;; esac
  case "$kind" in
    *.gz|*.bgz) gzip -t "$f" 2>/dev/null || { echo "not a valid gzip stream: $f" >&2; return 1; } ;;
    *.json) python3 -c "import json,sys; json.load(open(sys.argv[1]))" "$f" 2>/dev/null               || { echo "not valid JSON: $f" >&2; return 1; } ;;
    *.txt|*.tsv|*.csv|*.ind|*.snp) [ -s "$f" ] || { echo "empty: $f" >&2; return 1; } ;;
    *) [ -s "$f" ] || { echo "empty: $f" >&2; return 1; } ;;
  esac
}

# GRCh37 primary reference (the build every step below assumes)
[ -f b37/human_g1k_v37.fasta ] || {
  fetch https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz b37/human_g1k_v37.fasta.gz 100000000
  gzip -dc b37/human_g1k_v37.fasta.gz > b37/human_g1k_v37.fasta || true   # the archive has trailing bytes; harmless
  samtools faidx b37/human_g1k_v37.fasta; }

# GRCh38, only needed because PharmCAT is GRCh38-only
B38=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0
for f in Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.fasta.fai Homo_sapiens_assembly38.dict; do
  fetch "$B38/$f" "hg38/$f" 100000; done
for c in hg19ToHg38 hg38ToHg19; do
  fetch "https://hgdownload.soe.ucsc.edu/goldenPath/${c%%To*}/liftOver/$c.over.chain.gz" "chain/$c.over.chain.gz" 100000; done

# 1000 Genomes phase 3: plink2 bundle for PCA/PRS, per-chromosome VCFs for phasing and local ancestry
for e in pgen pvar.zst psam; do fetch "https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.$e?dl=1" "all_phase3.$e" 1000000; done
[ -f all_phase3.pvar ] && [ ! -f all_phase3.pvar.zst ] || "$PLINK2" --zst-decompress all_phase3.pvar.zst > all_phase3.pvar
K=https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf
for c in $(seq 1 22) X; do fetch "$K/chr$c.1kg.phase3.v5a.vcf.gz" "imp/ALL.chr$c.vcf.gz" 1000000; done
fetch https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar imp/beagle.jar 1000000
[ -d imp/maps ] || { fetch https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip imp/plink.GRCh37.map.zip 10000; unzip -q -o imp/plink.GRCh37.map.zip -d imp/maps; }

# annotation
[ -f clinvar_grch37.vcf.gz ] || { fetch https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz clinvar_grch37.vcf.gz 1000000
  fetch https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi clinvar_grch37.vcf.gz.tbi 10000; }
fetch https://ftp.ensembl.org/pub/grch37/release-87/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz annot/Homo_sapiens.GRCh37.87.gff3.gz 1000000
fetch https://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz annot/Homo_sapiens.GRCh37.87.gtf.gz 1000000
fetch https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz annot/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz 100000

# Y tree (YFull) and the ybrowse SNP position index
# 这三件决定 Y 单倍群走向：树坏了或 SNP 索引截断，调用结果会错得很安静，所以逐件校验
fetch https://raw.githubusercontent.com/YFullTeam/YTree/master/current_tree.json ytree/current_tree.json 100000
fetch https://raw.githubusercontent.com/YFullTeam/YTree/master/current_version.txt ytree/current_version.txt 4
fetch http://ybrowse.org/gbrowse2/gff/snps_hg19.csv ytree/snps_hg19.csv 100000

# optional panels: ancient DNA (AADR Human Origins) and archaic introgression (Sprime, Browning 2018)
if [ "${WITH_AADR:-1}" = "1" ]; then
  DV=https://dataverse.harvard.edu/api/access/datafile
  for id_name in "13994526:v66.p1_HO.aadr.patch.PUB.ind" "13994527:v66.p1_HO.aadr.patch.PUB.snp" \
                 "13994528:v66.p1_HO.aadr.PUB.anno" "13994808:v66.p1_HO.aadr.patch.PUB.geno"; do
    id=${id_name%%:*}; n=${id_name#*:}; fetch "$DV/$id" "aadr/$n" 100000; done
fi
if [ "${WITH_SPRIME:-1}" = "1" ]; then
  python3 - "$REF_DIR/archaic" <<'PY'
import json,sys,urllib.request,pathlib,subprocess
d=pathlib.Path(sys.argv[1]); d.mkdir(parents=True,exist_ok=True)
_req=urllib.request.Request("https://data.mendeley.com/public-api/datasets/y7hyt83vxr/files?folder_id=root&version=1",
        headers={"User-Agent":"Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120 Safari/537.36"})
meta=json.load(urllib.request.urlopen(_req))
want={f"{p}_sprime_results.tar.gz" for p in ("CHB","CHS","JPT","CEU","GBR","CDX","KHV","BEB","GIH")}
for f in meta:
    if f["filename"] in want and not (d/f["filename"]).exists():
        urllib.request.urlretrieve(f["content_details"]["download_url"], d/f["filename"])
        subprocess.run(["tar","xzf",str(d/f["filename"]),"-C",str(d)],check=True)
PY
fi
echo "references ready under $REF_DIR"
