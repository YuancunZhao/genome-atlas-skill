#!/usr/bin/env bash
# Run the whole pipeline in order. Any step that fails stops the run: do not work around a failure,
# fix it. Step 30 stops and waits for you to write work/report_text.yaml.
set -euo pipefail
cd "$(dirname "$0")"
source scripts/env.sh
S=scripts
# 每步的退出状态与产物关联（H6）：进入下一步时把上一步记为 ok，失败时由 trap ERR 记为 failed。
# 一次失败的运行必须在 run_info.json 里留下痕迹——否则下一轮会以为这一步跑过。
_STEP_NAME=""; _STEP_T0=$SECONDS
step(){
  if [ -n "$_STEP_NAME" ]; then
    python3 "$S/_run_info.py" --step "$_STEP_NAME" --status ok --seconds $((SECONDS-_STEP_T0)) >/dev/null 2>&1 || true
  fi
  echo -e "\n=== $* ==="; _STEP_NAME="$*"; _STEP_T0=$SECONDS
}
trap 'rc=$?; if [ -n "$_STEP_NAME" ]; then python3 "$S/_run_info.py" --step "$_STEP_NAME" --status failed --rc $rc >/dev/null 2>&1 || true; fi; echo "STEP FAILED (rc=$rc): $_STEP_NAME" >&2' ERR
# Ancestry steps follow the switches validated from config.yaml (read_options, §7 AN0). A step that is
# off writes a disabled manifest, so step 30 can tell "not configured" from "not run" instead of
# silently reusing an older result. A failing step still stops the run.
# 复审 AN0+AN5+H6：禁用记录原先写到 wgs/08、09、09b、16、17，而 30 检查的是 11_aadr、12_localanc、04_ancestry
# ——写的和读的不是同一处，于是"本次禁用 AADR"根本不会被看见，旧结果照旧进报告。现在写到**消费者实际
# 检查的目录**，文件名带步骤 id，避免同目录多步互相覆盖。
skipped(){ echo -e "\n=== skip $1: $2 ==="; python3 $S/ancestry_data.py --disabled "$1" \
  --out "$PROJ/wgs/$3/manifest.$1.json" --sample "$SAMPLE" --reason "$2"; }

# 运行事实（H6）：第一步之前记录 run id / 代码修订 / 有效参数 / 工具与参考版本。
# 事后回答不了「那次用的哪棵树、哪个参考、有没有未提交改动」是最常见的复现障碍。
python3 $S/_run_info.py || echo "warning: could not write run_info.json (continuing)"
step 00 prepare inputs and QC;                            bash  $S/00_qc.sh
step 01 normalise the VCF and build the callable mask;         bash  $S/01_normalize.sh
step 02 re-call X with the right ploidy and pile up MT;        bash  $S/02_recall_x_mt.sh
step 03 complete genotype set at the reference panel sites;    python3 $S/03_complete_set.py
step 04 ancestry PCA and projection;                           bash  $S/04_ancestry_pca.sh
step 04b nearest reference populations;                        python3 $S/04b_ancestry_summary.py
if [ "${REGIONAL_ENABLED:-0}" = "1" ] && [ -n "${AXIS:-}" ]; then
  step 04c per-chromosome axis index;                            python3 $S/04c_axis_index.py
  python3 $S/ancestry_data.py --clear-step 04c-per-chromosome-axis --dir "$PROJ/wgs/04_ancestry"
else
  skipped 04c-per-chromosome-axis regional_axis_not_configured 04_ancestry
fi
step 05 Y haplogroup on the YFull tree;                        python3 $S/05_y_haplogroup.py
step 06 mtDNA haplogroup and heteroplasmy;                     python3 $S/06_mtdna.py
step 06b mtDNA disease screen;                                 python3 $S/06b_mt_disease.py
step 07 annotate with ClinVar, consequences, frequencies;      bash  $S/07_annotate.sh
step 07b pathogenic and loss-of-function tables;               python3 $S/07b_clinvar_tables.py
step 07c gnomAD frequencies for the candidates;                python3 $S/07c_gnomad_lookup.py
# 医学与结构变异步骤（08a–08d、ROH）**不依赖祖源开关**：STR/SMN/SV/ROH 是报告医学结论的来源，
# 把它们放进 AADR 条件块意味着 AADR=0（默认）时这些结论整块不产出——曾如此，已移出。
step 08a main-contig alignment for SV calling;              bash  $S/08a_main_contigs.sh
step 08b repeat expansions;                                  bash  $S/08b_expansionhunter.sh
step 08c SMN1/SMN2 copy number;                            bash  $S/08c_smn.sh
step 08d structural variants via Delly, takes hours;          bash  $S/08d_delly.sh
step 09c runs of homozygosity;                                bash  $S/09c_roh.sh

if [ "${AADR_ENABLED:-0}" = "1" ]; then
  step 08 extract the ancient-DNA panel;                         python3 $S/08_aadr_extract.py
  step 09 ancient-DNA PCA and projection;                        bash  $S/09_aadr_pca.sh
  step 09b nearest present-day and ancient groups;               python3 $S/09b_aadr_summary.py
  # 复审 AN0/AN5/H6 生命周期：禁用轮写下的步骤级 disabled 记录，在对应步骤**本次成功**后必须撤下
  # （只撤自己的），否则 analysis_state 永远返回 disabled，重新启用等于没启用。
  python3 $S/ancestry_data.py --clear-step 08-aadr-extract      --dir "$PROJ/wgs/11_aadr"
  python3 $S/ancestry_data.py --clear-step 09-aadr-pca          --dir "$PROJ/wgs/11_aadr"
  python3 $S/ancestry_data.py --clear-step 09b-aadr-summary     --dir "$PROJ/wgs/11_aadr"
else
  skipped 08-aadr-extract aadr_not_configured 11_aadr
  skipped 09-aadr-pca aadr_not_configured 11_aadr
  skipped 09b-aadr-summary aadr_not_configured 11_aadr
fi
# Lineage history (AN4): consumes 05/06 results and the normalised metadata; runs before 30 so the
# report reads one file instead of re-deriving the paternal/maternal story from text.
step 09d lineage history and evidence;                        python3 $S/lineage_history.py --history panel/lineage_history.json --rows $WGS/11_aadr/summary.json --yard $WGS/03_haplo --out $WGS/03_haplo/lineage_history.json --sample "$SAMPLE"

step 10b CYP2D6 star alleles;                                bash  $S/10b_cyrius.sh
step 10 PharmCAT star alleles, CYP2D6, HLA;                    bash  $S/10_pharmcat.sh
step 11 extended pharmacogenomic markers;                      python3 $S/11_pgx_extra.py
step 11b HLA and KIR typing with T1K;                      bash  $S/11b_t1k.sh
step 12 polygenic scores;                                      python3 $S/12_prs.py
step 13 structural variants and copy number;                   python3 $S/13_sv_filter.py
step 14 read-backed phasing;                                   bash  $S/14_phase_reads.sh
step 14b phase summary and cis/trans questions;                python3 $S/14b_phase_summary.py
step 15 statistical phasing;                                   bash  $S/15_phase_statistical.sh
if [ "${LOCAL_ENABLED:-0}" = "1" ]; then
  step 16 local ancestry;                                        bash  $S/16_local_ancestry.sh
  # 校准是增强不是前置（AN3 复审）：16b 失败时原始 LA（17/17b）仍要独立交付。失败在这里降级——
  # run_info 记 failed，主流程继续。**不再**往 12_localanc 写 disabled 步骤记录：那份记录会被
  # analysis_state 当成整个 LA 目录的非 ok 状态，把本来可交付的原始结果一并吞掉（复审 AN0/AN5/H6）。
  # 16b 自己负责生命周期：开跑清旧记录与旧 calib.*，成功才写 state=ok 的凭证，17b 凭该记录读校准；
  # 失败时不留任何"本次校准成功"的痕迹，17b 把校准状态明确写成 not_run，报告不会把"没跑成校准"
  # 当成"跑了没写"，也不会再读上一轮的旧 calib.*。
  step 16b local-ancestry calibration
  if bash $S/16b_local_ancestry_calibration.sh; then
    :
  else
    _rc=$?
    python3 "$S/_run_info.py" --step "16b local-ancestry calibration" --status failed --rc "$_rc" >/dev/null 2>&1 || true
    _STEP_NAME=""
    echo "warning: 16b calibration failed (rc=$_rc); continuing -- raw local ancestry (17/17b) is still deliverable" >&2
  fi
  step 17 local-ancestry summary;                                python3 $S/17_local_ancestry_summary.py
  step 17b calibrated against held-out references;               python3 $S/17b_local_ancestry_calibrated.py
  python3 $S/ancestry_data.py --clear-step 16-local-ancestry    --dir "$PROJ/wgs/12_localanc"
  python3 $S/ancestry_data.py --clear-step 17-la-summary        --dir "$PROJ/wgs/12_localanc"
  python3 $S/ancestry_data.py --clear-step 17b-la-calibrated    --dir "$PROJ/wgs/12_localanc"
else
  skipped 16-local-ancestry local_ancestry_not_configured 12_localanc
  skipped 16b-la-calibration local_ancestry_not_configured 12_localanc
  skipped 17-la-summary local_ancestry_not_configured 12_localanc
  skipped 17b-la-calibrated local_ancestry_not_configured 12_localanc
fi
step 18 archaic introgressed segments;                         python3 $S/18_archaic.py
step 19 genes covered by those segments;                       python3 $S/19_archaic_genes.py
step 20 mutation spectrum;                                     python3 $S/20_mutation_spectrum.py
step 21 somatic signals in blood;                              python3 $S/21_somatic.py
step 22 telomere read fraction;                                bash  $S/22_telomere.sh
step 23 blood groups;                                          python3 $S/23_bloodgroups.py
step 24 KIR and HLA ligands;                                   python3 $S/24_kir.py
step 25 HLA disease and drug associations;                     python3 $S/25_hla_disease.py
step 26 behavioural polygenic scores;                          python3 $S/26_behaviour_prs.py
step 27 candidate behaviour genes;                             python3 $S/27_candidate_genes.py
step 28 f3 statistics;                                         python3 $S/28_f3_stats.py
step 30 assemble the report data;                              python3 $S/30_build_report_data.py

cat <<'MSG'

=== now write work/report_text.yaml ===
Copy example/report_text.yaml and rewrite every line for this sample. The scripts produce numbers;
the conclusions are yours to write, under the rules in SKILL.md. Then:

    python3 scripts/31_html_report.py

MSG
if [ -n "$_STEP_NAME" ]; then python3 "$S/_run_info.py" --step "$_STEP_NAME" --status ok --seconds $((SECONDS-_STEP_T0)) >/dev/null 2>&1 || true; fi
