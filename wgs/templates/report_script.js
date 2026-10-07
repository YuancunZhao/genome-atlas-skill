const D = __DATA__; const BLOODV3 = __BLOOD__; const UI = __UI__; const FIND = __FIND__; const PGX = __PGX__; const PRS_EN = __PRS_EN__;
let LANG = 'zh';
try { LANG = localStorage.getItem('dayu-lang') || (navigator.language && !navigator.language.startsWith('zh') ? 'en' : 'zh'); } catch(e) {}
const F = {};
const t = k => ((UI[k] || [k,k])[LANG === 'zh' ? 0 : 1]).replace(/\{(\w+)(:\+)?\}/g, (m,v,plus)=>{const x=F[v]; if(x===undefined) return m; return (plus && x>0 ? '+' : '') + x;});
const zh = () => LANG === 'zh';
const NAME = () => zh() ? (D.name_zh || 'Sample') : (D.name_en || 'Sample');
const TITLE = () => zh() ? (D.title_zh || D.name_zh || '基因组图鉴') : (D.title_en || D.name_en || 'Genome Atlas');
const NS='http://www.w3.org/2000/svg';
const el=(p,tag,a)=>{const n=document.createElementNS(NS,tag);for(const k in a)n.setAttribute(k,a[k]);p.appendChild(n);return n};
const txt=(p,a,s)=>{const n=el(p,'text',a);n.textContent=s;return n};
const tip=(n,s)=>{const q=document.createElementNS(NS,'title');q.textContent=s;n.appendChild(q)};
const V=name=>getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const C=()=>({ink:V('--ink'),muted:V('--muted'),lab:V('--lab'),faint:V('--faint'),grid:V('--grid'),track:V('--track'),data:V('--data'),data2:V('--data2'),fd:V('--faintdata'),hero:V('--hero'),bg:V('--bg')});
const fmt=n=>(n==null||Number.isNaN(n))?'—':n.toLocaleString('en-US');
const CH = {}; const NUM={'font-weight':800};
// 图注曾整体 toUpperCase()。全大写小字注是一种风格，但它会毁掉句子里的拉丁专名：中文图注里的
// SouthEA 会印成 SOUTHEA、Mb 印成 MB、call_rate 印成 CALL_RATE，看起来像编码而不是面板名。
// 保留字号/字距的注脚观感，去掉强制大写——信息正确优先于装饰。
const foot=(s,c,w,y,txt_)=>txt(s,{x:w/2,y,'text-anchor':'middle','font-size':7,'font-weight':600,fill:c.faint,'letter-spacing':'.12em',class:'fade',style:'animation-delay:.9s'},txt_);
// Clearing a figure with innerHTML='' is not portable: the containers are <svg> elements and
// Safari throws on writing innerHTML to an SVG element, which took the whole render down with it.
const clearEl=n=>{if(n.replaceChildren)n.replaceChildren();else while(n.firstChild)n.removeChild(n.firstChild)};
const reveal=(id,fn)=>{CH[id]=fn;const n=document.getElementById(id);n.style.cursor='pointer';n.addEventListener('click',()=>{clearEl(n);fn(n,C())})};
// HGNC convention: gene symbols are italic, protein names stay upright. The set comes from the data's
// own gene fields, never from a regex over the prose -- tokens such as PASS, CRAM, VCF, MGI, PLINK or
// the blood-group system names (ABO, Rh, FY ...) are not genes and must not be italicised.
const GENES=(()=>{const s=new Set();
  ['candidate','pgx','hla_disease','kir','sv_gene_dels','behaviour'].forEach(k=>(D[k]||[]).forEach(r=>{if(r&&r.gene)s.add(String(r.gene))}));
  (FIND||[]).forEach(f=>{if(f&&f.gene)s.add(String(f.gene))});   // findings live in the YAML payload, not in D
  (PGX||[]).forEach(r=>{if(r&&r.gene)s.add(String(r.gene))});    // as does the PGx table (CYP2D6, SLCO1B1 ...)
  ['ABO','Rh','FY','JK','MNS','KELL','LU','DI','DO','CO','LW','SC','IN','VEL','XG','YT','OK'].forEach(b=>s.delete(b));
  return s;})();
const gI=g=>'<i class="gname">'+g+'</i>';
// Wrap symbols in element text only (never inside tags), so marked names are not nested twice.
const gEm=h=>h.replace(/>([^<>]+)</g,(m,t)=>'>'+t.replace(/\b[A-Z][A-Z0-9]{2,11}\b/g,x=>GENES.has(x)?gI(x):x)+'<');
const K=D.kpi;
F.depth=K.depth_auto; F.x=K.depth_x; F.y=K.depth_y; F.mt=fmt(K.depth_mt); F.callable=K.callable_gb; F.pass=fmt(K.pass_records); F.snv=fmt(K.snv); F.indel=fmt(K.indel); F.titv=K.titv;
F.cmp=fmt(D.chip.compared); F.nocall=D.chip.nocall??'—'; F.nonp=D.chip.nocall_nopass??'—'; F.indel_n=D.chip.in_deletion??'—';
// AN5（§7）：KPI 与单倍群标量改从报告契约取（D.ancestry / D.lineages），而不是散落的派生键。
// 这里也修掉了缺 Y 时的无保护访问：D.ypath[D.ypath.length-1] 在没有 Y 结果时抛异常，
// 整个 renderAll 随之中断，20 张图一起变空白——一个缺失的模块不该带走整页。
const _AN=(D.ancestry&&D.ancestry.analyses)||[];
const _byScope=s=>_AN.find(a=>a.scope===s)||{};
const _byData=d=>_AN.find(a=>a.dataset===d)||{};
const _ranked=a=>analysisGroups(a).filter(g=>g.rank&&!g.small_group)
  .sort((x,y)=>(x.rank-y.rank)||(Number(x.distance_mean||0)-Number(y.distance_mean||0)));
const _sp=a=>Array.isArray(a.supported_path)?a.supported_path:[];
// 复审 P1（AN0/AN2/AN5）：百分位说明按**配置**的超群/亚群渲染（30 传 D.pop.superpop/subpops
// 与计数）。写死"东亚/汉"在 EUR 配置下数着 EUR 却自称东亚，是事实错误；未知代码原样显示。
const _SP_ZH={EAS:'东亚',EUR:'欧洲',AFR:'非洲',SAS:'南亚',AMR:'美洲'}, _SP_EN={EAS:'East Asians',EUR:'Europeans',AFR:'Africans',SAS:'South Asians',AMR:'Americans'};
const _spZh=()=>_SP_ZH[(D.pop||{}).superpop]||((D.pop||{}).superpop||'参考');
const _spEn=()=>_SP_EN[(D.pop||{}).superpop]||((D.pop||{}).superpop||'reference')+' individuals';
const _subPops=()=>{const ss=(D.pop||{}).subpops;return (Array.isArray(ss)&&ss.length)?ss:['CHB','CHS'];};
const _subZh=()=>_subPops().every(x=>x==='CHB'||x==='CHS')?'汉族':_subPops().join('/');
const _subEn=()=>_subPops().every(x=>x==='CHB'||x==='CHS')?'Han':_subPops().join('/');
F.ng=fmt(D.anc.n_global); F.ne=fmt(D.anc.n_eas);
F.knn=_ranked(_byScope('regional')).slice(0,3).map(g=>`${g.label} n=${g.n}`).join(' · ')||'—';
F.knng=_ranked(_byData('1000G')).slice(0,3).map(g=>`${g.label} n=${g.n}`).join(' · ')||'—';
F.near=(_ranked(_byScope('regional'))[0]||{}).label||(_ranked(_byData('1000G'))[0]||{}).label||'—';
F.nearg=_ranked(_byData('1000G')).slice(0,2).map(g=>`${g.label} ${(1*g.distance_mean).toFixed(4)}`).join(' ≈ ')||'—';
F.cvn=fmt(D.clinvar_total); F.cvd=D.clinvar_date; F.lof=D.lof.all; F.lofr=D.lof.rare; F.lofh=D.lof.rare_hom;
F.rho_i=D.prs_rho.imputed; F.rho_s=D.prs_rho.subset; F.sv=fmt(D.sv_total); F.rohn=D.roh_stats.n; F.rohmb=D.roh_stats.total_mb; F.rohmax=D.roh_stats.max_mb;
F.yterm=((D.lineages||{}).y||{}).reported_hg||D.y_terminal||"—"; F.yformed=(()=>{const p=_sp((D.lineages||{}).y||{})||[];return p.length?fmt(p[p.length-1].formed):"—";})(); F.mthg=((D.lineages||{}).mt||{}).reported_hg||((D.mt||{}).hg)||"—"; F.mtq=((D.lineages||{}).mt||{}).call_quality||((D.mt||{}).quality)||"—"; F.mtn=((D.mt||{}).found||[]).length; F.mtp=((D.mt||{}).private||[]).length;

/* ── AN6 父母系卡片：判定、分支时间线、折叠的背景层、发现记录与路线状态 ─────────────
   树的分叉不是地理迁移，formed/TMRCA 也不是迁移日期（§7）：时间线只画树上的时间，发现记录与
   路线分开列；没有来源路线时明确说"只有分布"，不为任何支系编故事。背景（祖先大支系）默认折叠。 */
const _esc=t=>String(t==null?'':t).replace(/[&<>"]/g,m=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[m]));
const renderLineageCards=()=>{
  const box=document.getElementById('linecard'); if(!box) return;
  const L=D.lineages||{}, Z=zh();
  const timeline=(path,width)=>{
    const pts=path.map(p=>({n:p.node,f:Number(p.formed),t:Number(p.tmrca)})).filter(p=>Number.isFinite(p.f));
    if(!pts.length) return '';
    const hi=Math.max(...pts.map(p=>p.f)), lo=0, X=v=>18+(width-36)*(1-Math.min(v,hi)/hi);
    const marks=pts.slice(-8).map(p=>`<line x1="${X(p.f).toFixed(1)}" y1="12" x2="${X(p.f).toFixed(1)}" y2="30" stroke="currentColor" stroke-width=".7" opacity=".5"/>`).join('');
    const ticks=[0,Math.round(hi/2),hi].map(v=>`<text x="${X(v).toFixed(1)}" y="44" font-size="7" text-anchor="middle" fill="currentColor" opacity=".6">${fmt(v)}</text>`).join('');
    return `<svg viewBox="0 0 ${width} 52" width="100%" height="52" role="img">
      <line x1="18" y1="30" x2="${width-18}" y2="30" stroke="currentColor" stroke-width="1" opacity=".45"/>${marks}${ticks}
      <text x="18" y="9" font-size="7" fill="currentColor" opacity=".6">${Z?'越靠左越古老（年）':'older (years) to the left'}</text></svg>`;
  };
  const card=(kind,lin,label_zh,label_en)=>{
    if(!lin) return '';
    const hg=_esc(lin.reported_hg||'—'), cons=lin.conservative_hg?_esc(lin.conservative_hg):'';
    const path=lin.supported_path||[], unc=(lin.uncertain_nodes||[]).map(u=>_esc(u.node||u));
    const hist=lin.history||{}, obs=hist.observations||[], routes=hist.routes||[];
    const state=_esc(lin.state||'ok'), why=_esc(lin.reason_code||'');
    const anc=path.filter(p=>!p.rank||true).slice(0,-4).map(p=>_esc(p.node)).filter(Boolean);
    const tail=path.slice(-4).map(p=>_esc(p.node));
    return `<div style="margin:10px 0 16px">
      <div style="font-family:var(--num);font-weight:700;font-size:13px">${Z?label_zh:label_en}: ${hg}${cons?` <span style="opacity:.75">（${Z?'保守回退':'conservative'} ${cons}）</span>`:''}</div>
      <div style="font-size:11px;opacity:.75;margin:2px 0 6px">${Z?'状态':'status'}: ${state}${why?' · '+why:''} ·
        ${Z?'树':'tree'} ${_esc(lin.tree_source||'')} ${_esc(lin.tree_version||'')} ·
        ${Z?'调用方式':'call'}: ${_esc(lin.call_quality||'')}</div>
      ${timeline(path, 420)}
      <div style="font-size:11px;margin-top:4px">${Z?'末端四级':'last four levels'}: ${tail.join(' → ')}${unc.length?` · ${Z?'支持位点不足':'weakly supported'}: ${unc.join('、')}`:''}</div>
      <details style="margin-top:6px"><summary style="cursor:pointer;font-size:11px;opacity:.8">${Z?'背景（更早的大支系）':'background (earlier major branches)'}</summary>
        <div style="font-size:11px;margin-top:4px;opacity:.8">${anc.join(' → ')||'—'}</div></details>
      <div style="font-size:11px;margin-top:6px">${Z?'已发表发现记录':'published records'}: ${obs.length}${obs.length?'':'（'+_esc(hist.history_reason_code||'')+'）'} ·
        ${Z?'迁移路线':'migration routes'}: ${routes.length}${routes.length?'':(Z?'（无有来源的路线，本卡只呈现分布）':' (no sourced route; distribution only)')}</div>
    </div>`;
  };
  box.innerHTML = card('y', L.y, '父系', "father's line") + card('mt', L.mt, '母系', "mother's line");
};

const renderAll=()=>{try{document.documentElement.setAttribute('lang',LANG==='zh'?'zh-Hans':'en');
  {const _bs=document.getElementById('bootstate');if(_bs)_bs.textContent='rendering…';}
  document.querySelectorAll('[data-i18n]').forEach(e=>{e.innerHTML=gEm(t(e.getAttribute('data-i18n')))});
  // The PCA legend had a literal `${NAME()}` in the static template -- JavaScript syntax that nothing
  // substitutes, so it reached the page as text. Fill it here with the other runtime strings, which
  // also makes it follow the language toggle.
  {const _lm=document.getElementById('lbl_me');if(_lm)_lm.textContent=NAME();}
  document.getElementById('langbtn').textContent=t('toggle');
  document.title = TITLE();
  window.__renderErrors=[];
  for(const id in CH){const n=document.getElementById(id);if(!n)continue;clearEl(n);
    try{CH[id](n,C())}catch(e){window.__renderErrors.push(id);console.error('[render] '+id,e);
      const err=document.createElementNS('http://www.w3.org/2000/svg','text');
      err.setAttribute('x',8);err.setAttribute('y',16);err.setAttribute('font-size',10);err.setAttribute('fill','#b3261e');
      err.textContent=id+': '+(e&&e.message?e.message:e);n.appendChild(err)}}
  if(window.__renderErrors.length)console.warn('[render] failed figures:',window.__renderErrors.join(', '));
  renderKpi(); renderFindings(); renderPgx(); renderHla(); renderMisc(); renderSections(); renderLineageCards(); renderF3(); if(typeof deepRender==='function') deepRender();{const _bd=document.getElementById('bootstate');if(_bd){const _nf=(window.__renderErrors||[]).length,_nt=Object.keys(CH).length;_bd.textContent='rendered · figures '+(_nt-_nf)+'/'+_nt+(_nf?' (failed: '+window.__renderErrors.join(', ')+')':'');_bd.style.color=_nf?'#ffb4a2':'#9ae6b4';}}}catch(e){console.error('[renderAll]',e);{const _be=document.getElementById('bootstate');if(_be){_be.textContent='RENDER FAILED: '+(e&&e.message?e.message:e);_be.style.color='#ffb4a2';}}try{document.body.insertAdjacentHTML('afterbegin','<div style="margin:14px auto;max-width:900px;padding:10px 14px;border:1px solid #b3261e;border-radius:8px;color:#b3261e;font:13px/1.5 system-ui">render error: '+(e&&e.message?e.message:e)+'</div>')}catch(_){}}};
document.getElementById('langbtn').addEventListener('click',()=>{LANG=zh()?'en':'zh';try{localStorage.setItem('dayu-lang',LANG)}catch(e){};renderAll()});
if(window.matchMedia){window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change',renderAll)}
const prsName=n=>zh()?(PRS_EN[n]||n):n;

function renderKpi(){
  document.getElementById('kpi_depth').innerHTML=K.depth_auto+'<small>×</small>';
  document.getElementById('kpi_var').innerHTML=(K.pass_records/1e6).toFixed(2)+'<small>M</small>';
  document.getElementById('kpi_call').innerHTML=K.callable_gb+'<small>Gb</small>';
  document.getElementById('kpi_titv').textContent=K.titv.toFixed(2);
}

/* ── 00 per-chromosome depth + PASS variants ── */
reveal('chromdepth',(s,c)=>{
  const rows=D.chrom, x0=50, top=14, rh=12.2, W1=300, W2=260;
  txt(s,{x:x0,y:6,'font-size':7,'font-weight':600,fill:c.faint,'letter-spacing':'.08em'},zh()?'平均深度 0–40×':'MEAN DEPTH 0–40×');
  txt(s,{x:x0+W1+120,y:6,'font-size':7,'font-weight':600,fill:c.faint,'letter-spacing':'.08em'},zh()?'PASS 变异数（每 Mb）':'PASS VARIANTS PER MB');
  const maxd=Math.max(...rows.filter(r=>r.chrom!=='MT').map(r=>r.depth)), maxr=Math.max(...rows.filter(r=>r.len>0&&r.chrom!=='MT').map(r=>r.n_pass/r.len*1e6));
  rows.forEach((r,i)=>{const y=top+i*rh+6, sex=['X','Y','MT'].includes(r.chrom);
    txt(s,{x:x0-8,y:y+3,'text-anchor':'end','font-size':8,'font-weight':700,fill:c.ink,'letter-spacing':'.06em'},r.chrom);
    el(s,'line',{x1:x0,y1:y,x2:x0+W1,y2:y,stroke:c.grid,'stroke-width':.6});
    const d=r.chrom==='MT'?40:r.depth, bw=W1*Math.min(d,40)/40;
    const b=el(s,'rect',{x:x0,y:y-3.6,width:bw,height:7.2,rx:1.5,fill:sex?c.data2:c.data,class:'fade',style:`animation-delay:${i*.03}s`}); tip(b,`chr${r.chrom}: ${r.depth}×`);
    if(r.chrom==='MT') el(s,'rect',{x:x0+W1-4,y:y-3.6,width:2,height:7.2,fill:c.bg});
    txt(s,{x:x0+bw+5,y:y+3,'font-size':8,'font-weight':800,fill:c.ink},r.depth+'×');
    const px=x0+W1+120, ratio=r.len>0?r.n_pass/r.len*1e6:0, bw2=W2*Math.min(ratio,maxr)/maxr;
    el(s,'line',{x1:px,y1:y,x2:px+W2,y2:y,stroke:c.grid,'stroke-width':.6});
    const b2=el(s,'rect',{x:px,y:y-3.6,width:r.chrom==='MT'?0:bw2,height:7.2,rx:1.5,fill:sex?c.data2:c.data,class:'fade',style:`animation-delay:${.2+i*.03}s`}); tip(b2, r.chrom==='MT' ? (zh()?'chrMT：PASS 调用集不含 MT，变异见母系分析':'chrMT: the PASS call set has no MT records; see the MT panel') : `chr${r.chrom}: ${fmt(r.n_pass)} PASS variants`);
    txt(s,{x:px+(r.chrom==='MT'?0:bw2)+5,y:y+3,'font-size':7.5,'font-weight':800,fill:c.muted},r.chrom==='MT'?(zh()?'见母系分析':'see MT panel'):fmt(r.n_pass));
  });
  [0,10,20,30,40].forEach(v=>{const x=x0+W1*v/40,yb=top+rows.length*rh+8;el(s,'line',{x1:x,y1:yb,x2:x,y2:yb+4,stroke:c.faint,'stroke-width':.8});txt(s,{x,y:yb+13,'text-anchor':'middle','font-size':7,fill:c.faint},v+'×')});
  foot(s,c,900,348,zh()?`X 约为常染色体一半、Y 更低：男性 · 线粒体 ${fmt((D.chrom.find(r=>r.chrom==='MT')||{}).depth||0)}× 超出刻度 · MT 的变异由独立 pileup 分析给出（PASS 调用集不含 MT，故此处不计）`:`X about half of autosomes, Y lower: male · mitochondria ${fmt((D.chrom.find(r=>r.chrom==='MT')||{}).depth||0)}× off scale · MT variants come from a separate pileup analysis (the PASS call set has no MT records, so none are counted here)`);
});


/* ── 01 Y chain (YFull path walked by 05 from the tree root; any haplogroup) ── */
reveal('ychain',(s,c)=>{
  // Only the last six steps are drawn, with the upstream levels collapsed into a dotted stub. The
  // full 24-level walk put an 831-site root node next to 1-site tips, which crushed every label,
  // and those early levels carry no conclusion: the sample's story is the resolution limit at the
  // tip. The bundled example report draws it the same way (six circles plus an upstream text line).
  const FULL=D.ypath, nF=FULL.length, TAIL=Math.min(6,nF);
  const P=FULL.slice(-TAIL), n=P.length, up=FULL.slice(0,nF-TAIL);
  const y=100, base=170, x1=490, x0=40;
  el(s,'line',{x1:base,y1:y,x2:x1,y2:y,stroke:c.ink,'stroke-width':1,class:'draw'});
  if(up.length){
    el(s,'circle',{cx:52,cy:y,r:5,fill:c.data,opacity:.35});
    el(s,'line',{x1:58,y1:y,x2:base-10,y2:y,stroke:c.faint,'stroke-width':1.4,'stroke-dasharray':'1.5 4'});
    txt(s,{x:52,y:y-13,'text-anchor':'middle','font-size':11,'font-weight':700,fill:c.faint},'⋯');
    txt(s,{x:52,y:y+21,'text-anchor':'middle','font-size':7.5,'font-weight':700,fill:c.muted},up.length+(zh()?' 级上游':' upstream'));
    txt(s,{x:52,y:y+32,'text-anchor':'middle','font-size':6.5,fill:c.faint},up[0].snp);
  }
  // Radius is sqrt-normalised into a fixed range. Area still scales with the derived-site count,
  // which is what the caption promises.
  const _maxDer=Math.max(1,...P.map(p=>p.der));
  const _rad=der=>7+20*Math.sqrt(der/_maxDer);
  P.forEach((p,i)=>{const x=base+(x1-base)*i/(n-1), r=_rad(p.der), last=i===n-1;
    const cir=el(s,'circle',{cx:x,cy:y,r,fill:last?c.hero:c.data,stroke:c.bg,'stroke-width':1.5,class:'pop'}); cir.style.animationDelay=(i*110)+'ms';
    // der>0: counted derived sites. der=0 with anc>0: genuinely ancestral. Both zero: the tree lists no
    // hg19-mapped site for this branch, which is not the same as a coverage failure.
    const _st = p.der>0 ? (p.der+(zh()?' 位点':' SNP'))
             : (p.anc>0 ? (zh()?`0 衍/${p.anc} 祖`:`0 der/${p.anc} anc`)
             : (p.na>0 ? (zh()?`覆盖不足 ${p.na}`:`low cov ${p.na}`) : (zh()?'无 hg19 位点':'no hg19 site')));
    tip(cir,`${p.snp} · der=${p.der} anc=${p.anc} n/a=${p.na} · formed ${fmt(p.formed)} ybp`);
    txt(s,{x,y:y-r-10,'text-anchor':'middle','font-size':last?11.5:10,'font-weight':last?800:700,fill:c.ink,stroke:c.bg,'stroke-width':2.6,style:'paint-order:stroke',...NUM},p.snp.replace('O-',''));
    txt(s,{x,y:y+r+14,'text-anchor':'middle','font-size':7.5,'font-weight':600,fill:c.muted,stroke:c.bg,'stroke-width':2.6,style:'paint-order:stroke'},_st);
    txt(s,{x,y:y+r+26,'text-anchor':'middle','font-size':7,fill:c.faint,stroke:c.bg,'stroke-width':2.6,style:'paint-order:stroke'},'~'+fmt(p.formed)+(zh()?' 年前':' ybp'));
  });
  txt(s,{x:x0,y:y+80,'font-size':8,'font-weight':600,fill:c.muted},(D.ypath&&D.ypath.length)?(zh()?`YFull 树路径（共 ${D.ypath.length} 级）：${D.ypath.slice(0,8).map(p=>p.snp).join(' → ')} → … → ${D.ypath[D.ypath.length-1].snp}`:`YFull path (${D.ypath.length} levels): ${D.ypath.slice(0,8).map(p=>p.snp).join(' → ')} → … → ${D.ypath[D.ypath.length-1].snp}`):(zh()?'Y 路径数据缺失':'Y path data unavailable'));
  txt(s,{x:x0,y:y+96,'font-size':8,'font-weight':600,fill:c.muted},zh()?'逐级从 YFull 树下行，每级取衍生态支持最多的支系（明细见 y_terminal_snps.tsv）':'Walked down the YFull tree, taking the best derived-state-supported branch at each step');
  txt(s,{x:x0,y:y+130,'font-size':24,'font-weight':700,fill:c.ink,...NUM},D.y_terminal);
  txt(s,{x:x0+150,y:y+130,'font-size':9,'font-weight':600,fill:c.muted},'YFull '+(D.versions.yfull||'—')+(zh()?' · 叶节点':' · leaf branch'));
  const _ysD=(D.y_snps||[]).map(r=>r.depth).filter(d=>d>0);
  const _ysState=(D.y_snps||[]).reduce((a,r)=>{a[r.state]=(a[r.state]||0)+1;return a;},{});
  const _ysMin=_ysD.length?Math.min(..._ysD):0, _ysMax=_ysD.length?Math.max(..._ysD):0;
  const _ysMed=_ysD.length?[..._ysD].sort((a,b)=>a-b)[Math.floor(_ysD.length/2)]:0;
  txt(s,{x:x0,y:y+150,'font-size':8,'font-weight':600,fill:c.muted},_ysD.length?(zh()?`${_ysD.length} 个末端定义位点（${Object.entries(_ysState).map(([k,v])=>k+' '+v).join(' · ')}）· 实测深度 ${_ysMin}–${_ysMax}×（中位 ${_ysMed}×）`:`${_ysD.length} terminal defining sites (${Object.entries(_ysState).map(([k,v])=>k+' '+v).join(' · ')}) · observed depth ${_ysMin}-${_ysMax}x (median ${_ysMed}x)`):(zh()?'末端定义位点数据缺失':'terminal defining-site data unavailable'));
  txt(s,{x:x0,y:y+163,'font-size':8,'font-weight':600,fill:c.faint},zh()?'质控门：MAPQ≥30 · BQ≥25 · DP≥5 · 读段一致性≥90%（深度不足的位点单列为 nocov，不并入衍生态计数）':'quality gates: MAPQ>=30 · BQ>=25 · DP>=5 · read concordance>=90% (under-covered sites are listed as nocov, never counted as derived)');
  foot(s,c,520,296,zh()?'图上只画末端 6 级（更早 18 级见下方路径）· 圆面积 = 该级衍生态位点数 · 最深色 = 终端 · 每级位点均通过上述质控门':'only the last six steps are drawn (the earlier 18 are listed below) · circle area = derived sites at that step · darkest = terminal · every step passes the gates above');
});

/* ── 01 mt strip ── */
reveal('mtstrip',(s,c)=>{
  // 复审 P1（AN0/AN5）：无 mt 是交付形态（未做 mt 调用/只有 Y 的样本），画明确说明而不是
  // 在 D.mt.found 处带着整页崩溃。
  if(!D.mt||!D.mt.hg){
    txt(s,{x:200,y:140,'text-anchor':'middle','font-size':12,'font-weight':700,fill:c.muted},zh()?'线粒体结果不可用':'mt result unavailable');
    txt(s,{x:200,y:162,'text-anchor':'middle','font-size':8.5,fill:c.faint},
        zh()?'未找到单倍群判定产物（haplogrep3.txt）——母系结论按不可用交付，不用常染色体相似性补位'
            :'no haplogroup call was produced (haplogrep3.txt absent); the maternal line is delivered as unavailable, never backfilled from autosomal similarity');
    return;
  }
  const L=16569, x0=24, x1=376, y=110;
  el(s,'line',{x1:x0,y1:y,x2:x1,y2:y,stroke:c.ink,'stroke-width':1});
  [0,4000,8000,12000,16569].forEach(p=>{const x=x0+(x1-x0)*p/L;el(s,'line',{x1:x,y1:y+4,x2:x,y2:y+8,stroke:c.faint,'stroke-width':.8});txt(s,{x,y:y+18,'text-anchor':'middle','font-size':7,fill:c.faint},p===16569?'16,569':fmt(p))});
  const hvs=x0+(x1-x0)*16024/L; el(s,'rect',{x:hvs,y:y-30,width:x1-hvs,height:26,fill:c.fd,opacity:.5}); txt(s,{x:x1,y:y-34,'text-anchor':'end','font-size':7,'font-weight':600,fill:c.muted},'HVS-I');
  el(s,'rect',{x:x0,y:y-30,width:(x1-x0)*576/L,height:26,fill:c.fd,opacity:.5}); txt(s,{x:x0,y:y-34,'font-size':7,'font-weight':600,fill:c.muted},'D-loop');
  const all=D.mt.found.map(v=>[v,false]).concat(D.mt.private.map(v=>[v,true])).sort((a,b)=>parseFloat(a[0])-parseFloat(b[0]));
  all.forEach(([v,isP],i)=>{const pos=parseFloat(v), x=x0+(x1-x0)*pos/L, ln=i%2===0?-22:22;
    el(s,'line',{x1:x,y1:y,x2:x,y2:y+ln*.55,stroke:c.faint,'stroke-width':.7});
    const d=el(s,'circle',{cx:x,cy:y+ln*.55,r:3.2,fill:isP?c.bg:c.data,stroke:isP?c.ink:c.data,'stroke-width':1,class:'pop'}); d.style.animationDelay=(i*20)+'ms'; tip(d,v+(isP?(zh()?' · 私有':' · private'):(zh()?' · 定义位点':' · defining')));
  });
  const _h0=(D.mt_het||[])[0]; if(_h0){const hx=x0+(x1-x0)*_h0.pos/L; el(s,'path',{d:`M${hx-5} ${y+34} L${hx} ${y+27} L${hx+5} ${y+34} Z`,fill:c.hero}); txt(s,{x:hx,y:y+45,'text-anchor':'middle','font-size':7,'font-weight':700,fill:c.ink},`${_h0.pos}${_h0.alt||''} ${(_h0.af*100).toFixed(0)}%`)}
  txt(s,{x:x0,y:y+80,'font-size':22,'font-weight':700,fill:c.ink,...NUM},D.mt.hg);
  const _mtF=D.mt.found.length;
  const _mtTag = _mtF ? (zh()?`${_mtF} 个定义位点命中`:`${_mtF} defining sites`)
                      : (zh()?'HSD 内无变异：PASS VCF 不含 MT 调用，单倍群仅据 pileup 覆盖判断':`no variants in the HSD: the PASS VCF has no MT calls, so the call rests on pileup coverage alone`);
  txt(s,{x:x0+80,y:y+80,'font-size':10,'font-weight':600,fill:c.muted},'PhyloTree '+D.versions.phylotree+' · '+_mtTag+(zh()?' · 质量 ':' · quality ')+D.mt.quality);
  const _hetTxt=(D.mt_het||[]).map(h=>`m.${h.pos}${h.alt||''} ${(h.af*100).toFixed(0)}%`).join('、');
  const _hetNote=(D.mt_het||[]).length && (D.mt_het||[]).every(h=>[310,3107].includes(h.pos));
  const lines=zh()?['未命中的支系定义位点：'+D.mt.notfound.join('、'),'私有变异（'+D.mt.private.length+'）：'+D.mt.private.join(' '),'异质性：'+(_hetTxt||'未检出')+(_hetNote?'（310/3107 是已知参考序列伪影位点）':'')]:
    ['Defining sites not found: '+D.mt.notfound.join(', '),'Private ('+D.mt.private.length+'): '+D.mt.private.join(' '),'Heteroplasmy: '+(_hetTxt||'none detected')+(_hetNote?' (m.310/3107 are known reference artefacts)':'')];
  lines.forEach((tx,i)=>txt(s,{x:x0,y:y+104+i*14,'font-size':7.8,'font-weight':500,fill:c.muted},tx));
  foot(s,c,400,296,zh()?`实心 = 支系定义位点 · 空心 = 私有变异 · 三角 = 异质性位点 · ${fmt((D.chrom.find(r=>r.chrom==='MT')||{}).depth||0)}× 深度`:`filled = branch-defining · hollow = private · triangle = heteroplasmic · ${fmt((D.chrom.find(r=>r.chrom==='MT')||{}).depth||0)}× depth`);
});

/* ── 02 scatter ── */
function scatter(s,c,pts,me,opts){
  const {x0=44,y0=14,W=360,H=290,xl='PC1',yl='PC2',color,labels}=opts;
  const xs=pts.map(p=>p.x).concat(me.x), ys=pts.map(p=>p.y).concat(me.y);
  let xmin=Math.min(...xs),xmax=Math.max(...xs),ymin=Math.min(...ys),ymax=Math.max(...ys);
  const dx=(xmax-xmin)*.06,dy=(ymax-ymin)*.06; xmin-=dx;xmax+=dx;ymin-=dy;ymax+=dy;
  const X=v=>x0+W*(v-xmin)/(xmax-xmin), Y=v=>y0+H*(1-(v-ymin)/(ymax-ymin));
  el(s,'rect',{x:x0,y:y0,width:W,height:H,fill:'none',stroke:c.grid,'stroke-width':.6});
  txt(s,{x:x0+W/2,y:y0+H+22,'text-anchor':'middle','font-size':8,'font-weight':700,fill:c.muted,'letter-spacing':'.1em'},xl);
  txt(s,{x:x0-30,y:y0+H/2,'text-anchor':'middle','font-size':8,'font-weight':700,fill:c.muted,'letter-spacing':'.1em',transform:`rotate(-90 ${x0-30} ${y0+H/2})`},yl);
  pts.forEach(p=>{const d=el(s,'circle',{cx:X(p.x).toFixed(1),cy:Y(p.y).toFixed(1),r:2.1,fill:color(p),opacity:.85});tip(d,p.t)});
  if(labels){const groups={};pts.forEach(p=>{(groups[p.g]=groups[p.g]||[]).push(p)});
    Object.entries(groups).forEach(([g,arr])=>{const mx=arr.reduce((a,p)=>a+p.x,0)/arr.length,my=arr.reduce((a,p)=>a+p.y,0)/arr.length;
      txt(s,{x:X(mx),y:Y(my)-6,'text-anchor':'middle','font-size':7.5,'font-weight':700,fill:c.ink,stroke:c.bg,'stroke-width':2.5,'paint-order':'stroke','letter-spacing':'.06em'},g.toUpperCase())})}
  const mx=X(me.x),my=Y(me.y);
  el(s,'polygon',{points:`${mx},${my-7.8} ${mx+7.8},${my} ${mx},${my+7.8} ${mx-7.8},${my}`,fill:c.hero,stroke:c.bg,'stroke-width':1.5,class:'pop'});
  txt(s,{x:mx+9,y:my+3,'font-size':9,'font-weight':800,fill:c.ink,stroke:c.bg,'stroke-width':3,'paint-order':'stroke'},NAME());
  if(opts.foot) foot(s,c,420,336,opts.foot);
}
reveal('pcaglobal',(s,c)=>{const pts=D.pca_global.pts.map(p=>({g:p[0],x:p[2],y:p[3],t:`${p[1]} (${p[0]})`}));
  scatter(s,c,pts,{x:D.pca_global.me[0],y:D.pca_global.me[1]},{color:p=>p.g==='EAS'?c.data:c.fd,labels:true,foot:zh()?'一点 = 一个人 · 深蓝 = 东亚 · 菱形 = '+NAME():'one dot = one person · dark blue = East Asian · diamond = '+NAME()})});
reveal('pcaeas',(s,c)=>{
  // 复审 P1（AN0/AN2/AN5）：区域视图按状态交付。未配置 ref_superpop 是合法的 global-only
  // 形态、配置了但产物缺失是待重跑——两者都画明确说明，不造空图也不拿另一份参考冒充。
  if(!D.pca_eas){
    const r=D.pca_eas_reason;
    txt(s,{x:210,y:130,'text-anchor':'middle','font-size':12,'font-weight':700,fill:c.muted},
        r==='missing_products'?(zh()?'区域参考视图不可用':'regional reference view unavailable')
                            :(zh()?'区域参考视图未配置':'regional reference view not configured'));
    txt(s,{x:210,y:152,'text-anchor':'middle','font-size':8.5,fill:c.faint},
        r==='missing_products'
          ?(zh()?'ref_superpop 已配置，但 04 的 regional.* 产物缺失——重跑 04_ancestry_pca.sh 后可恢复':'ref_superpop is configured but 04\'s regional.* products are missing; re-run 04_ancestry_pca.sh to restore the view')
          :(zh()?'ref_superpop 未给出：本报告按全球参考交付，区域视图不是缺件':'no ref_superpop was given: this report is delivered on the global reference; the regional view is not a missing piece'));
    return;
  }
  const shade={CHB:c.hero,CHS:c.data,JPT:c.data2,CDX:c.fd,KHV:c.fd};const pts=D.pca_eas.pts.map(p=>({g:p[0],x:p[1],y:p[2],t:p[0]}));
  scatter(s,c,pts,{x:D.pca_eas.me[0],y:D.pca_eas.me[1]},{color:p=>shade[p.g]||c.fd,labels:true,foot:zh()?'一点 = 一个人 · 明度 = 人群 · 菱形 = '+NAME():'one dot = one person · shade = population · diamond = '+NAME()})});

/* ── 03 findings ── */
function renderFindings(){
  document.getElementById('findings').innerHTML=FIND.map(f=>`<div class="finding"><div><div class="g">${gI(f.gene)}<small>${zh()?f.gt_zh:f.gt_en}</small></div><span class="stamp${f.hot?' hot':''}">${zh()?f.tag_zh:f.tag_en}</span></div><div><div class="b">${zh()?f.zh:f.en}</div>${f.act_zh?`<div class="act">${zh()?f.act_zh:f.act_en}</div>`:''}</div></div>`).join('');
}
/* ── 03 ClinVar bars ── */
reveal('cvbars',(s,c)=>{
  const order=['Benign','Likely_benign','Uncertain_significance','Conflicting_classifications_of_pathogenicity','drug_response','association','risk_factor','other','not_provided','protective','Pathogenic/Likely_pathogenic','Pathogenic','Likely_pathogenic'];
  const rows=order.filter(k=>D.clinvar[k]).map(k=>[k,D.clinvar[k]]); const x0=250,x1=880,top=10,rh=13, max=rows[0][1];
  const nm={Benign:['良性','Benign'],Likely_benign:['可能良性','Likely benign'],Uncertain_significance:['意义未明','Uncertain'],Conflicting_classifications_of_pathogenicity:['分类冲突','Conflicting'],drug_response:['药物反应','Drug response'],association:['关联','Association'],risk_factor:['风险因素','Risk factor'],other:['其它','Other'],not_provided:['未提供','Not provided'],protective:['保护性','Protective'],'Pathogenic/Likely_pathogenic':['致病/可能致病','P/LP'],Pathogenic:['致病','Pathogenic'],Likely_pathogenic:['可能致病','Likely pathogenic']};
  rows.forEach(([k,v],i)=>{const y=top+i*rh+5, w=Math.max(1.5,(x1-x0)*Math.sqrt(v)/Math.sqrt(max)), pat=/athogenic/.test(k)&&!/Conflicting/.test(k);
    txt(s,{x:x0-8,y:y+3,'text-anchor':'end','font-size':8,'font-weight':pat?800:600,fill:pat?c.ink:c.muted},(nm[k]||[k,k])[zh()?0:1]);
    const b=el(s,'rect',{x:x0,y:y-3.5,width:w,height:7,rx:1.5,fill:pat?c.hero:(/Conflicting|Uncertain/.test(k)?c.data:c.data2),class:'fade',style:`animation-delay:${i*.04}s`}); tip(b,`${k}: ${fmt(v)}`);
    txt(s,{x:x0+w+5,y:y+3,'font-size':8,'font-weight':800,fill:c.ink},fmt(v));
  });
  foot(s,c,900,186,zh()?'条长按平方根压缩 · 最深色 = 致病/可能致病 · 只统计 PASS 变异':'bar length on a square-root scale · darkest = pathogenic / likely pathogenic · PASS variants only');
});

/* ── 04 PGx ── */
function renderPgx(){
  document.getElementById('pgx').innerHTML=PGX.map(r=>`<tr><td class="gene">${gI(r.gene)}</td><td class="g">${r.dip}</td><td>${zh()?r.ph_zh:r.ph_en}</td><td>${zh()?r.zh:r.en}</td><td class="rs">${r.src}</td></tr>`).join('');
}
function renderHla(){
  const order=['HLA-A','HLA-B','HLA-C','HLA-DRB1','HLA-DQA1','HLA-DQB1','HLA-DPA1','HLA-DPB1'];
  const hot=new Set();
  document.getElementById('hla').innerHTML=order.map(g=>{const v=D.hla[g]; if(!v) return ''; const f=a=>hot.has(a)?`<em>${a.replace(/:\d+$/,'')}</em>`:a.replace(/(\*\d+:\d+).*/,'$1'); return `<div><div class="k">${g}</div><div class="v">${f(v[0])} / ${f(v[1])}</div></div>`}).join('');
}

/* ── 05 PRS three versions ── */
reveal('prs',(s,c)=>{
  const rows=[...D.prs].sort((a,b)=>b.pct_EAS-a.pct_EAS), x0=250, x1=850, top=16, rh=21, X=v=>x0+(x1-x0)*v/100;
  // The panel grows with the number of scores, but the canvas was a fixed 700 tall while 47 rows
  // need 16+47*21 = 1003: the caption landed in the middle of the data and the last rows spilled
  // out of the canvas into the report text below. Size the canvas from the content instead.
  const H=top+rows.length*rh+44; s.setAttribute('viewBox','0 0 900 '+H);
  el(s,'rect',{x:X(25),y:top-8,width:X(75)-X(25),height:rows.length*rh+4,fill:c.track,opacity:.6});
  [0,25,50,75,100].forEach(v=>{el(s,'line',{x1:X(v),y1:top-8,x2:X(v),y2:top+rows.length*rh-2,stroke:v===50?c.faint:c.grid,'stroke-width':v===50?.9:.5});txt(s,{x:X(v),y:top+rows.length*rh+12,'text-anchor':'middle','font-size':7.5,fill:c.faint},v)});
  rows.forEach((r,i)=>{const y=top+i*rh+6, hi=r.pct_EAS>=90||r.pct_EAS<=10;
    txt(s,{x:x0-12,y:y+3,'text-anchor':'end','font-size':9,'font-weight':hi?800:600,fill:hi?c.ink:c.lab},prsName(r.trait));
    el(s,'line',{x1:X(r.pct_Han),y1:y,x2:X(r.pct_EAS),y2:y,stroke:c.data,'stroke-width':1,opacity:.5}); const q=el(s,'circle',{cx:X(r.pct_Han),cy:y,r:3.6,fill:c.bg,stroke:c.data,'stroke-width':1.3});tip(q,(zh()?'汉族 208 人中 ':'among 208 Han ')+r.pct_Han)
    const mx=X(r.pct_EAS); const d=el(s,'polygon',{points:`${mx},${y-6.4} ${mx+6.4},${y} ${mx},${y+6.4} ${mx-6.4},${y}`,fill:c.hero,stroke:c.bg,'stroke-width':1.2,class:'pop'}); d.style.animationDelay=(i*25)+'ms'; tip(d,`${r.trait}: WGS ${r.pct_EAS} (Han ${r.pct_Han}) · coverage ${r.coverage_pct}% · z ${r.z_vs_EAS}`);
    txt(s,{x:x1+10,y:y+3,'font-size':8.5,'font-weight':800,fill:c.ink},r.pct_EAS.toFixed(0));
    txt(s,{x:x1+34,y:y+3,'font-size':7.5,fill:c.faint},r.coverage_pct+'%');
  });
  txt(s,{x:x1+10,y:top-10,'font-size':7,'font-weight':600,fill:c.faint,'letter-spacing':'.08em'},zh()?'百分位 · 覆盖':'PCT · COVERAGE');
  foot(s,c,900,H-8,zh()?`菱形 = 在 ${(D.pop||{}).n_super||'—'} 名${_spZh()}人中的百分位 · 空心圆 = 在 ${(D.pop||{}).n_sub||'—'} 名${_subZh()}中 · 灰带 = 中间一半的人`:`diamond = percentile among ${(D.pop||{}).n_super||'—'} ${_spEn()} · hollow = among ${(D.pop||{}).n_sub||'—'} ${_subEn()} · band = middle half`);
});

/* ── 06 SV panel (Delly): the slot used to hold a copy-number depth chart that no producer fills,
      while the calls that do exist (3,619 SV, 10 whole-gene deletion events) went unshown. ── */
reveal('cnv',(s,c)=>{
  const rows=D.sv_by_chrom||[], gd=D.sv_gene_dels||[], st=D.sv_gene_dels_stats||{}, cnt=D.sv_counts||{};
  if(!rows.length){
    txt(s,{x:12,y:22,'font-size':10,fill:c.faint},zh()?'本版没有结构变异产物（Delly 未产出）。':'No structural-variant calls in this build (Delly produced nothing).');return}
  const x0=52, x1=270, top=32, rh=11.2, maxc=Math.max(1,...rows.map(r=>r.DEL+r.DUP+r.INV)), X=v=>x0+(x1-x0)*v/maxc;
  txt(s,{x:x0-4,y:16,'font-size':8,'font-weight':600,fill:c.faint,'letter-spacing':'.06em'},zh()?'每条染色体上的 SV（堆叠）':'SV CALLS PER CHROMOSOME');
  rows.forEach((r,i)=>{const y=top+i*rh; let acc=0;
    txt(s,{x:x0-6,y:y+3,'text-anchor':'end','font-size':7.5,fill:c.faint},r.chrom);
    [['DEL',c.data],['DUP',c.data2],['INV',c.faint]].forEach(([k,col])=>{const v=r[k]||0; if(!v)return;
      const bx=X(acc), bw=Math.max(X(acc+v)-bx,.8);
      const b=el(s,'rect',{x:bx,y:y-3.2,width:bw,height:6.4,rx:1,fill:col,class:'fade',style:`animation-delay:${i*.02}s`});
      tip(b,`chr${r.chrom}: ${k} ${v}`); acc+=v});
    txt(s,{x:X(acc)+5,y:y+3,'font-size':7,fill:c.faint},acc)});
  const ly=top+rows.length*rh+16;
  [['DEL',c.data],['DUP',c.data2],['INV',c.faint]].forEach(([k,col],i)=>{
    el(s,'rect',{x:x0+i*92,y:ly-5.5,width:9,height:7,rx:1,fill:col});
    txt(s,{x:x0+i*92+13,y:ly+1,'font-size':7.5,fill:c.muted},k+' '+fmt(cnt[k]||0))});
  const rx=300, rtop=42, rrh=13, BW=38;
  txt(s,{x:rx,y:16,'font-size':8,'font-weight':600,fill:c.faint,'letter-spacing':'.06em'},zh()?'整段缺失的基因':'WHOLE-GENE DELETIONS');
  txt(s,{x:rx,y:27,'font-size':7.5,fill:c.muted},(st.events||0)+(zh()?' 个事件 / ':' events / ')+(st.genes||0)+(zh()?' 个基因 · 深度核验':' genes · depth-checked'));
  gd.slice(0,14).forEach((r,i)=>{const y=rtop+i*rrh, f=(typeof r.frac==='number')?r.frac:null;
    txt(s,{x:rx,y:y+3,'font-size':8,'font-weight':700,fill:c.ink},r.gene);
    txt(s,{x:rx+110,y:y+3,'font-size':7.5,fill:c.faint},'chr'+r.chrom);
    el(s,'line',{x1:rx+128,y1:y,x2:rx+128+BW,y2:y,stroke:c.grid,'stroke-width':3.4});
    if(f!==null) el(s,'line',{x1:rx+128,y1:y,x2:rx+128+BW*Math.max(.04,Math.min(f,1)),y2:y,stroke:c.hero,'stroke-width':3.4});
    txt(s,{x:rx+128+BW+6,y:y+3,'font-size':7.5,'font-weight':700,fill:c.ink},f===null?'—':f.toFixed(2)+'×')});
  foot(s,c,520,324,zh()?`Delly 2.6.0 读段对判定 + 深度核验 · 共 ${fmt(D.sv_total||0)} 个 SV（大小中位 ${D.sv_size_median||'—'} bp）· 右侧条 = 该基因上读段深度 / 全基因组平均，越低越像真缺失`:`Delly 2.6.0 read-pair calls + depth checks · ${fmt(D.sv_total||0)} SVs (median size ${D.sv_size_median||'—'} bp) · right-hand bar = read depth over the gene / genome average, lower = more deletion-like`);
});
/* ── 06 STR panel ── */
reveal('str',(s,c)=>{
  const rows=D.str.filter(r=>r.thr).sort((a,b)=>b.max/b.thr-a.max/a.thr).slice(0,22), x0=78, x1=372, top=14, rh=13.6, X=v=>x0+(x1-x0)*Math.min(v,1.1)/1.1;
  el(s,'line',{x1:X(1),y1:top-6,x2:X(1),y2:top+rows.length*rh,stroke:c.ink,'stroke-width':1}); txt(s,{x:X(1),y:top-9,'text-anchor':'middle','font-size':7,'font-weight':700,fill:c.ink},zh()?'自身阈值':'own threshold');
  rows.forEach((r,i)=>{const y=top+i*rh+6, w=X(r.max/r.thr)-x0;
    txt(s,{x:x0-6,y:y+3,'text-anchor':'end','font-size':8,'font-weight':700,fill:c.ink},r.locus);
    el(s,'line',{x1:x0,y1:y,x2:X(1),y2:y,stroke:c.grid,'stroke-width':.5});
    const b=el(s,'rect',{x:x0,y:y-3.5,width:Math.max(w,1),height:7,rx:1.5,fill:(r.max>=r.thr)?'#d2644a':c.data,class:'fade',style:`animation-delay:${i*.03}s`}); tip(b,`${r.locus} (${r.unit}): ${r.gt} repeats${r.ci?` (CI ${r.ci})`:''} · threshold ≥${r.thr}`);
    txt(s,{x:x0+w+4,y:y+3,'font-size':7.5,'font-weight':800,fill:c.ink},r.gt);
  });
  // The caption states what the data shows (closest locus, its ratio, any locus at/above the
  // recorded threshold, and how many assayed loci carry no local rule) -- it never asserts
  // "all below threshold" unconditionally.
  const all=D.str||[], ruled=all.filter(r=>r.thr), over=ruled.filter(r=>r.max>=r.thr);
  const worst=ruled.reduce((a,b)=>(!a||b.max/b.thr>a.max/a.thr)?b:a,null);
  const norule=all.filter(r=>!r.thr).length;
  foot(s,c,400,326,zh()
    ?`条长 = 最长等位基因/该位点自身阈值（逐位点不同）${worst?` · 最接近：${worst.locus} ${(worst.max/worst.thr).toFixed(2)}×`:''}${over.length?` · ⚠ ${over.length} 个位点已达阈值`:' · 本样本无位点达阈值'}${norule?` · 另有 ${norule} 个已测位点无本地规则、未纳入本图`:''} · 置信区间见悬停 · ExpansionHunter 5`
    :`bar = longest allele / recorded threshold${worst?` · closest: ${worst.locus} ${(worst.max/worst.thr).toFixed(2)}×`:''}${over.length?` · ⚠ ${over.length} locus/loci at or above threshold`:' · no locus at threshold in this sample'}${norule?` · ${norule} assayed loci have no local rule and are not shown`:''} · confidence intervals on hover · ExpansionHunter 5`);
});

/* ── 07 ROH ideogram ── */
reveal('roh',(s,c)=>{
  const chr=[...Array(22)].map((_,i)=>String(i+1)), x0=40, top=12, rh=12.4, maxL=D.chrlen['1'], W=830;
  chr.forEach((ch,i)=>{const y=top+i*rh+5, L=D.chrlen[ch], w=W*L/maxL;
    txt(s,{x:x0-8,y:y+3,'text-anchor':'end','font-size':8,'font-weight':700,fill:c.ink},ch);
    el(s,'rect',{x:x0,y:y-2.5,width:w,height:5,rx:2.5,fill:c.track});
    const cen=D.cen[ch]; if(cen) el(s,'line',{x1:x0+W*cen*1e6/maxL,y1:y-5,x2:x0+W*cen*1e6/maxL,y2:y+5,stroke:c.faint,'stroke-width':.9});
    D.roh.filter(r=>r.chrom===ch).forEach(r=>{const rx=x0+W*r.start/maxL, rw=Math.max(2,W*(r.end-r.start)/maxL); const b=el(s,'rect',{x:rx,y:y-3.2,width:rw,height:6.4,rx:1,fill:r.mb>=1.5?c.hero:c.data,class:'pop'}); tip(b,`chr${ch}:${fmt(r.start)}-${fmt(r.end)} · ${r.mb} Mb · q ${r.q}`)});
  });
  foot(s,c,900,296,zh()?'深色段 = ≥1.5 Mb 的纯合区 · 浅色 = 1–1.5 Mb · 小刻度 = 着丝粒 · 跨着丝粒/端粒的伪影段已剔除':'dark = homozygous run ≥1.5 Mb · light = 1–1.5 Mb · tick = centromere · segments spanning centromeres/telomeres removed');
});
// Consumes D.sections (step 30's status table): every section not "ok" is listed with its
// machine reason code and bilingual reason, so an empty figure is explained instead of
// reading as a negative result. All-ok reports render nothing here.
function renderSections(){
  const el=document.getElementById('secstatus'); if(!el) return;
  const bad=(D.sections||[]).filter(s=>s&&s.status!=='ok');
  if(!bad.length){el.innerHTML='';return;}
  const L=s=>zh()?(s.name_zh||s.id):(s.name_en||s.id);
  const R=s=>zh()?(s.reason_zh||s.code||''):(s.reason_en||s.code||'');
  el.innerHTML=(zh()
    ?'<b>数据状态：</b>以下模块本次不可用，对应的空图或空表<b>不是阴性结果</b>——<br>'
    :'<b>Data status:</b> the sections below were unavailable; their empty figures are <b>not negative results</b> —<br>')
    +bad.map(s=>`· ${L(s)} — ${R(s)}${s.evidence?` <span style="opacity:.55">(${s.evidence})</span>`:''}`).join('<br>');
}
// SMNCopyNumberCaller's isCarrier means "SMN1 copy number == 1" -- it is not an exclusion of
// every carrier configuration. Two SMN1 copies can still sit in cis (2+0), and variants the
// caller does not model are untouched, so a false flag never prints "not a carrier".
const smnCell=()=>{
  if(!(D.smn&&D.smn.SMN1!=null)) return zh()?'未评估（无 SMN 拷贝数产物）':'not assessed (no SMN copy-number output)';
  const b=zh()?`${D.smn.SMN1} / ${D.smn.SMN2} 份`:`${D.smn.SMN1} / ${D.smn.SMN2} copies`;
  if(D.smn.SMN1===1) return b+(zh()?' · SMN1 单拷贝（携带者）':' · single SMN1 copy (carrier)');
  return b+(zh()?` · SMN1 ${D.smn.SMN1} 拷贝：未检出单拷贝型，但顺式 2+0 构型与本方法未建模的变异不能排除，残余携带风险仍存在`:` · SMN1 ${D.smn.SMN1} copies: no single-copy type detected, but a cis 2+0 arrangement and variants this caller does not model are not excluded, so residual carrier risk remains`);
};
function renderMisc(){  const _roh5=zh()?((D.roh_stats.n_gt5?`${D.roh_stats.n_gt5} 段 >5 Mb`:'无 >5 Mb')):((D.roh_stats.n_gt5?`${D.roh_stats.n_gt5} >5 Mb`:'none >5 Mb'));  const items=zh()?[['SMN1 / SMN2',smnCell()],['线粒体异质性',(D.mt_het||[]).map(h=>`m.${h.pos}${h.alt||''} ${(h.af*100).toFixed(0)}%`).join('、')||'未检出'],['结构变异总数',`${fmt(D.sv_total)}（DEL ${fmt(D.sv_counts.DEL)} · DUP ${D.sv_counts.DUP} · INV ${D.sv_counts.INV}）`],['罕见功能丧失变异',`${D.lof.rare} 个 · 纯合 ${D.lof.rare_hom} 个（重复区状态见变异表）`],['X 染色体',(D.par_het!=null?`PAR 重调后 ${fmt(D.par_het)} 个杂合位点`:'PAR 重调数据缺失')],['近亲检查',`ROH ≥1 Mb ${D.roh_stats.n} 段 · 合计 ${D.roh_stats.total_mb} Mb · ${_roh5}`]]:
    [['SMN1 / SMN2',smnCell()],['mt heteroplasmy',(D.mt_het||[]).map(h=>`m.${h.pos}${h.alt||''} ${(h.af*100).toFixed(0)}%`).join(', ')||'none detected'],['Structural variants',`${fmt(D.sv_total)} (DEL ${fmt(D.sv_counts.DEL)} · DUP ${D.sv_counts.DUP} · INV ${D.sv_counts.INV})`],['Rare loss-of-function',`${D.lof.rare} · ${D.lof.rare_hom} homozygous (repeat status in the variant table)`],['X chromosome',(D.par_het!=null?`${fmt(D.par_het)} heterozygous PAR sites after re-calling`:'PAR re-call data unavailable')],['Relatedness',`${D.roh_stats.n} ROH ≥1 Mb · ${D.roh_stats.total_mb} Mb total · ${_roh5}`]];
  document.getElementById('misc').innerHTML=items.map(([k,v])=>`<div><div class="k">${k}</div><div class="v" style="font-size:12px">${v}</div></div>`).join('');
}

/* ══════════ deep-dive additions (phase 6) ══════════ */
const AC=()=>({n:V('--ancn'),s:V('--ancs'),a:V('--arch')});
const CHR22=[...Array(22)].map((_,i)=>String(i+1));
/* AN6（§7）：校准面板、人数与染色体一律来自结构化结果，不再写死 CHB/CHS 与 1/2/6/22。
   配置里换一组来源面板或换几条染色体，这里会跟着变；拿不到就以空数组降级（图上不出现假值）。 */
const _LOC=(D.ancestry&&D.ancestry.local)||{};
const _CAL=Array.isArray(_LOC.calibration)?_LOC.calibration:[];
const _SRC=((_LOC.panels||[]).find(p=>p.role==='source')||{});
// 参与比较的面板来自校准结果本身（配置的 la_labels），不是 FLARE 的列顺序。
// 复审 AN3/AN5：没有校准时，标签回退到**原始 LA** 已知的来源面板（17 写入 panels），不能变 '—'。
const _SRC_PANELS=(_LOC.panels||[]).filter(p=>p.role==='source').map(p=>p.id||p.label).filter(Boolean);
const _PANELS=laPanelLabels(_LOC.calibration_panels,_LOC.panels);
const _SRCLAB=_PANELS[0]||_SRC.label||_SRC.id||'';
const CALCH=(Array.isArray(_LOC.calibration_chroms)&&_LOC.calibration_chroms.length)
  ? _LOC.calibration_chroms.map(String) : [];
/* 复审 AN3/AN5：KPI 与校准图必须是同一个统计量。dayuCal 取 17b 在共同染色体集合上写回的
   posterior 目标值（target_north），不再从 17 的 per_chrom 片段跨度比例折算——片段比例 0.9 而
   posterior 目标 0.6/参照均值 0.5/SD 0.1 时，旧代码显示 90%/+4 SD，同口径的正解是 60%/+1 SD。
   校准行没有 target_north（17b 未运行）时 dayuCal 为 null：KPI 显示 '—'，不用别的口径顶上。 */
function calibTarget(cal){
  const r=(cal||[]).find(x=>x&&x.target_north!=null&&isFinite(Number(x.target_north)));
  return r?Number(r.target_north):null;
}
function laSdScore(t,mean,sd){
  if(t==null||mean==null||!(Number(sd)>0)) return null;
  if(!isFinite(t)||!isFinite(mean)||!isFinite(sd)) return null;
  return (t-mean)/sd;
}
/* 复审 AN3/AN5：空 calibration_panels 不能把标签变 '—'——原始 LA 的 panels（17 写入）已经
   声明了来源面板；校准在场时仍以校准面板（配置的 la_labels）为准。 */
function laPanelLabels(calPans, rawPans){
  const cal=(calPans||[]).filter(Boolean);
  if(cal.length) return cal;
  return (rawPans||[]).filter(p=>p&&p.role==='source').map(p=>p.id||p.label).filter(Boolean);
}
const dayuCal=calibTarget(_CAL);
/* AN6（§7）：局部祖源的面板名、对照面板与"哪一列是来源"都来自结构化结果，不再写死
   NorthEA/SouthEA/European/SouthAsian。拿不到就留 '—'——不编一个看起来合理的百分比。
   _PANELS 已在上方带原始 LA 回退。 */
const _lg=D.la_global||{};
const _P0=_PANELS[0]||'', _P1=_PANELS[1]||'';
const _CTRL=((_LOC.panels||[]).filter(p=>p.role==='control').map(p=>p.id||p.label)).filter(Boolean);
F.north=(_P0&&_lg[_P0]!=null)?(_lg[_P0]*100).toFixed(1):'—';
F.south=(_P1&&_lg[_P1]!=null)?(_lg[_P1]*100).toFixed(1):'—';
F.la_a=_P0||'—'; F.la_b=_P1||'—';
F.noise=(_CTRL.length&&_CTRL.every(k=>_lg[k]!=null))?(_CTRL.reduce((a,k)=>a+_lg[k],0)*100).toFixed(1):'—';
F.chb=_CAL[0]&&_CAL[0].north_mean!=null?(_CAL[0].north_mean*100).toFixed(0):'—';
F.chs=_CAL[1]&&_CAL[1].north_mean!=null?(_CAL[1].north_mean*100).toFixed(0):'—';
F.dayucal=dayuCal!=null?(dayuCal*100).toFixed(1):'—';
const _sds=laSdScore(dayuCal,_CAL[0]&&_CAL[0].north_mean,_CAL[0]&&_CAL[0].north_sd);
F.sdchb=_sds==null?'—':((_sds>0?'+':'')+ _sds.toFixed(1));
F.archmb=D.archaic_summary.span_mb; F.neamb=D.archaic_summary.neanderthal_mb; F.denmb=D.archaic_summary.denisovan_mb;
F.archn=D.archaic_summary.merged; F.archpct=(D.archaic_summary.span_mb/2875*100).toFixed(1);
F.phhet=fmt(D.phase.het); F.phased=fmt(D.phase.phased); F.phpct=D.phase.het>0?(D.phase.phased/D.phase.het*100).toFixed(1):'—';
F.n50=(D.phase.n50_kb/1000).toFixed(1); F.blocks=fmt(D.phase.blocks);
F.mtcn=D.somatic.mtDNA_copies_per_cell; F.ydr=D.somatic.Y_depth_ratio; F.telk7=fmt(D.telomere.k7);
F.teltot=fmt((D.telomere&&D.telomere.total_reads)||0); F.archhom=((D.archaic_summary||{}).homozygous)||0;
F.phmax=(((D.phase||{}).max_mb)||0).toFixed(1); F.chipalt=fmt(((D.chip_hotspots||{}).alt_reads)||0);
F.segmax=Math.max(0,...((D.la_segments||[]).filter(g=>_P1&&g.anc===_P1).map(g=>g.mb))).toFixed(1);
F.snvn=fmt(D.spectrum.reduce((a,b)=>a+b.n,0));
F.cpg=(D.spectrum.filter(x=>/\[C>T\]G/.test(x.ctx)).reduce((a,b)=>a+b.frac,0)*100).toFixed(1);

const ANC_ZH={"China_Henan_Jiaozuoniecunsite_LBA_IA":"河南 焦作聂村 晚商–铁器时代","China_Shandong_Chengziya_Yueshi":"山东 城子崖 岳石文化","China_Henan_Pingliangtaisite_LN":"河南 平粮台 龙山文化","China_Baligang_BA_EasternZhou":"河南 八里岗 东周","China_Jinan_LiuJiaZhuang_Shang":"山东 济南刘家庄 商","China_Jining_YinJiaCheng_Longshan":"山东 济宁尹家城 龙山","China_Shandong_Chengziya_Longshan":"山东 城子崖 龙山","China_InnerMongolia_Erdaojingzi_LN":"内蒙古 二道井子 新石器晚期","China_Henan_Haojiatai_LN":"河南 郝家台 龙山","China_Baligang_LN_Longshan":"河南 八里岗 龙山","China_Shandong_Dinggong_LN":"山东 丁公 龙山","China_Henan_Wadiansite_LN":"河南 瓦店 龙山","China_Baligang_LN_Shijiahe":"河南 八里岗 石家河","China_Qingdao_BeiQian_Dawenkou":"山东 青岛北阡 大汶口","China_Baligang_LN_Yangshao":"河南 八里岗 仰韶","China_Shanxi_Shengedaliang_LN":"陕西 神圪垯梁 龙山","China_LBA_EIA":"中国 晚青铜–早铁器","China_Qinghai_Dacaozisite_IA":"青海 大草子 铁器时代","China_Qinghai_Lajiasite_LN":"青海 喇家 龙山","Japan_KofunPeriod":"日本 古坟时代",
  // 补全：报告的"所属群体"列与散点标签会对所有出现过的群体取中文名，缺条目就会漏出英文 ID
  "China_MLBA":"中国 中晚期青铜（含云南白羊村 M13）","China_MBA":"中国 中期青铜","China_LBA":"中国 晚期青铜",
  "China_IA":"中国 铁器时代","China_IA_Hellenistic":"中国 铁器时代（希腊化期）","China_LN_Xiaoheyan":"中国 小河沿 新石器晚期",
  "Taiwan_EN":"台湾 新石器早期","Taiwan_IA":"台湾 铁器时代",
  "China_AmurRiverBasin_LatePaleolithic":"黑龙江流域 旧石器晚期","China_AmurRiverBasin_Mesolithic":"黑龙江流域 中石器","China_AmurRiverBasin_N":"黑龙江流域 新石器",
  "China_Qinghai_Zongri":"青海 宗日",
  "China_Tibet_Agangrong":"西藏 阿岗绒","China_Tibet_Chaxiutang":"西藏 茶秀塘","China_Tibet_Piyangjiweng":"西藏 皮央吉翁",
  "China_Tibet_Gebusailu_IA":"西藏 格布赛鲁 铁器时代","China_Tibet_Butaxiongqu":"西藏 布塔雄曲","China_Tibet_Latuotanggu":"西藏 拉托塘果",
  "China_Tibet_Pulanduowa_IA":"西藏 普兰多瓦 铁器时代","China_Tibet_Longsangquduo":"西藏 陇桑曲多","China_Tibet_Qulongsazha_IA":"西藏 曲龙萨扎 铁器时代",
  "China_Tibet_Ounie":"西藏 欧聂","China_Tibet_Nudagang":"西藏 努达岗","China_Tibet_Sangdalongguo_IA":"西藏 桑达隆果 铁器时代",
  "Mongolia_XiongnuPeriod":"蒙古 匈奴时期","Mongolia_East_N":"蒙古 东部 新石器","Mongolia_N":"蒙古 新石器",
  "Japan_Nagabaka_2800BP":"日本 长坂 2800 BP","Japan_Chiba_HG_Jomon":"日本 千叶 绳文（狩猎采集）",
  "Japan_Honshu_EarlyJomon":"日本 本州 绳文早期","Japan_Honshu_MidLateJomon":"日本 本州 绳文中晚期"};

/* ── chromosome painting ── */
reveal('painting',(s,c)=>{
  const a=AC(), x0=34, rows=11, rh=44, maxL=D.chrlen['1'], COLW=390;
  const seg={}; D.la_segments.forEach(g=>{(seg[g.chrom]=seg[g.chrom]||[]).push(g)});
  const _miss=new Set((D.la_missing||[]).map(String));
  CHR22.forEach((ch,i)=>{
    const col=i<rows?0:1, row=i%rows, bx=col===0?x0:x0+470, y=18+row*rh, w=COLW*D.chrlen[ch]/maxL;
    txt(s,{x:bx-6,y:y+11,'text-anchor':'end','font-size':10,'font-weight':700,fill:c.ink},ch);
    [0,1].forEach(h=>{const yy=y+h*11;
      el(s,'rect',{x:bx,y:yy,width:w,height:9,rx:1.5,fill:a.n,opacity:.5});
      (seg[ch]||[]).filter(g=>g.hap===h+1).forEach(g=>{
        const gx=bx+w*g.start/D.chrlen[ch], gw=Math.max(1.2,w*(g.end-g.start)/D.chrlen[ch]);
        // 片段色必须与图例色块同源（lg_n 色块 = 实心 --ancn）：面板 A 实心主色（数据层画成
        // 页面底色 c.bg 等于不可见）、面板 B 强调色、对照面板中性色（它们只是对照）
        const r=el(s,'rect',{x:gx,y:yy,width:gw,height:9,rx:1,
                             fill:(_P1&&g.anc===_P1)?a.s:(g.anc===_P0?a.n:c.fd)});
        tip(r,`chr${ch}:${fmt(g.start)}-${fmt(g.end)} · ${g.anc}`)})});
    const cen=D.cen[ch]; if(cen) el(s,'line',{x1:bx+w*cen*1e6/D.chrlen[ch],y1:y-2,x2:bx+w*cen*1e6/D.chrlen[ch],y2:y+22,stroke:c.bg,'stroke-width':1.6});
    // A chromosome with no segments has no inference at all (FLARE left a 0-byte file and step 17
    // dropped it in silence). Say so, instead of leaving a band that reads like the genuinely untyped
    // acrocentric short arms.
    if(_miss.has(ch)){
      el(s,'rect',{x:bx,y:y,width:w,height:20,rx:1.5,fill:'#000',opacity:.07});
      txt(s,{x:bx+w/2,y:y+13,'text-anchor':'middle','font-size':8,'font-weight':700,fill:c.muted},zh()?'无推断：该染色体未产出（数据缺失）':'no inference: nothing produced for this chromosome');
    } else {
      const so=D.la_per_chrom.find(r=>String(r.chrom)===ch);
      if(so&&_P1) txt(s,{x:bx+w+7,y:y+13,'font-size':9,'font-weight':700,fill:c.muted},(so[_P1]*100).toFixed(0)+'%');
    }
  });
  // 面板名来自结果（calibration_panels），不是写死的"北方/南方"——换一对来源面板，这句话跟着换
  foot(s,c,900,496,zh()?`每条染色体两行 = 父母各给的一条单倍型 · 只画 ≥0.5 Mb 的「${F.la_b}」片段 · 右侧 = 该染色体的「${F.la_b}」比例 · 灰条 = 无推断（数据缺失，不是 0）· 短臂留白 = 参考面板在该区无标记`
                       :`two rows per chromosome = the two parental haplotypes · only ≥0.5 Mb segments of "${F.la_b}" · right = that chromosome's "${F.la_b}" share · grey = no inference (missing data, not zero) · blank short arms = no panel markers there`);
});

/* ── calibration ── */
reveal('calib',(s,c)=>{
  const a=AC(), x0=170,x1=770,y0=34,rh=32;
  if(!_CAL.length){txt(s,{x:12,y:22,'font-size':10,fill:c.faint},
    zh()?'本版没有留出校准结果（未运行或来源面板不足）':'No holdout calibration in this build (not run, or too few reference individuals)');return}
  // 每个校准群体一行：名字、人数、均值与标准差都来自结果本身。样本行画的是 17b 的 posterior
  // 目标值（与各组均值同口径）；没有 target_north 就不画样本行，不拿片段跨度比例冒充。
  const rows=(dayuCal!=null?[{lab:NAME(),n:dayuCal,sd:0,me:true}]:[]).concat(
    _CAL.map(x=>({lab:`${x.population} · n=${x.n}`,n:x.north_mean,sd:x.north_sd||0})));
  const X=v=>x0+(x1-x0)*(v-0.5)/0.5;
  [0.5,0.6,0.7,0.8,0.9,1].forEach(v=>{el(s,'line',{x1:X(v),y1:y0-12,x2:X(v),y2:y0+rows.length*rh-10,stroke:c.grid,'stroke-width':.7});
    txt(s,{x:X(v),y:y0+rows.length*rh+6,'text-anchor':'middle','font-size':8,fill:c.faint},(v*100)+'%')});
  txt(s,{x:x0,y:14,'font-size':8,'font-weight':600,fill:c.faint,'letter-spacing':'.08em'},
    (zh()?`${CALCH.join('、')} 号染色体上的「${_SRCLAB||'来源面板'}」成分`
        :`"${_SRCLAB||'source panel'}" share on chromosome(s) ${CALCH.join(', ')}`));
  rows.forEach((r,i)=>{const y=y0+i*rh;
    txt(s,{x:x0-12,y:y+4,'text-anchor':'end','font-size':10.5,'font-weight':r.me?800:600,fill:r.me?c.ink:c.lab},r.lab);
    if(r.sd>0){el(s,'line',{x1:X(r.n-r.sd),y1:y,x2:X(r.n+r.sd),y2:y,stroke:c.faint,'stroke-width':1.4});
      [-1,1].forEach(k=>el(s,'line',{x1:X(r.n+k*r.sd),y1:y-4,x2:X(r.n+k*r.sd),y2:y+4,stroke:c.faint,'stroke-width':1.4}))}
    const d=el(s,'circle',{cx:X(r.n),cy:y,r:r.me?6:4.5,fill:r.me?a.n:c.bg,stroke:r.me?c.bg:c.lab,'stroke-width':r.me?2:1.5,class:'pop'});
    tip(d,`${(r.n*100).toFixed(1)}%${r.sd?` ± ${(r.sd*100).toFixed(1)}`:''}`);
    txt(s,{x:X(r.n+(r.sd||0))+15,y:y+4,'font-size':10.5,'font-weight':800,fill:c.ink},(r.n*100).toFixed(1)+'%')});
  foot(s,c,900,166,zh()?'横线 = 组内 ±1 标准差 · 参考个体打分时已从参考面板移出':'whisker = ±1 sd within the group · reference individuals held out of the panel when scored');
});

/* ── Human Origins PCA with ancient projection ── */
reveal('hopca',(s,c)=>{
  const a=AC(), x0=44,y0=16,W=440,H=390;
  const pts=D.ho_modern, anc=D.ho_ancient_pts, me=D.ho_me;
  const xs=pts.map(p=>p[1]).concat(anc.map(p=>p.pc1),[me[0]]), ys=pts.map(p=>p[2]).concat(anc.map(p=>p.pc2),[me[1]]);
  let xmin=Math.min(...xs),xmax=Math.max(...xs),ymin=Math.min(...ys),ymax=Math.max(...ys);
  const dx=(xmax-xmin)*.05,dy=(ymax-ymin)*.05; xmin-=dx;xmax+=dx;ymin-=dy;ymax+=dy;
  const X=v=>x0+W*(v-xmin)/(xmax-xmin), Y=v=>y0+H*(1-(v-ymin)/(ymax-ymin));
  el(s,'rect',{x:x0,y:y0,width:W,height:H,fill:'none',stroke:c.grid,'stroke-width':.6});
  // AN6：环上的颜色按"来源面板 A / 来源面板 B / 其它"分配，名字来自结果而不是固定成分含义
  const COL={}; COL[_P0||'A']=a.n; COL[_P1||'B']=a.s; COL.none=c.fd;
  pts.forEach(p=>{el(s,'circle',{cx:X(p[1]).toFixed(1),cy:Y(p[2]).toFixed(1),r:2,fill:COL[p[0]]||c.fd,opacity:COL[p[0]]?.6:.34})});
  anc.forEach(p=>{const d=el(s,'rect',{x:X(p.pc1)-2.4,y:Y(p.pc2)-2.4,width:4.8,height:4.8,fill:'none',stroke:c.ink,'stroke-width':1,opacity:.75});
    tip(d,`${p.label} · ${fmt(p.date)} BP`)});
  // D.ho_prov carries one row per sample, so labelling every row drew 67 labels and 67 leader lines on
  // top of each other. Group by province instead, the way a clustering figure does it: one colour per
  // province, one label at the group centroid, and a colour key in the corner.
  const PROV_ZH={Shandong:'山东',Henan:'河南',Fujian:'福建',Guangdong:'广东',Sichuan:'四川',
                 Hubei:'湖北',Jiangsu:'江苏',Shanxi:'陕西',Zhejiang:'浙江',Shanghai:'上海',Chongqing:'重庆'};
  const PROV_C={Shandong:'#2f6fb0',Henan:'#c0504d',Fujian:'#4f9d69',Guangdong:'#8e6bbf',
                Sichuan:'#d08c2a',Hubei:'#3aa6a6',Jiangsu:'#b0558e',Shanxi:'#7f8f2a',
                Zhejiang:'#c76f9c',Shanghai:'#7a8b99',Chongqing:'#a0522d'};
  const _pg={}; (D.ho_prov||[]).forEach(p=>{(_pg[p.label]=_pg[p.label]||[]).push(p)});
  const _pcs=Object.entries(_pg).map(([k,arr])=>({label:k,n:arr.length,
      x:arr.reduce((s,p)=>s+p.pc1,0)/arr.length, y:arr.reduce((s,p)=>s+p.pc2,0)/arr.length}))
    .sort((p,q)=>q.n-p.n);
  _pcs.forEach(p=>{const col=PROV_C[p.label]||c.ink, px=X(p.x), py=Y(p.y);
    el(s,'circle',{cx:px,cy:py,r:3.2,fill:col,stroke:c.bg,'stroke-width':1.3,class:'pop'});
    txt(s,{x:px,y:py-9,'text-anchor':'middle','font-size':9,'font-weight':700,fill:col,stroke:c.bg,'stroke-width':2.6,'paint-order':'stroke'},
        zh()?(PROV_ZH[p.label]||p.label):p.label);
    txt(s,{x:px,y:py+14,'text-anchor':'middle','font-size':6.5,'font-weight':600,fill:c.muted,stroke:c.bg,'stroke-width':2.2,'paint-order':'stroke'},'n='+p.n)});
  {let kx=x0+8, ky=y0+8;
   _pcs.slice(0,6).forEach((p,i)=>{const col=PROV_C[p.label]||c.ink;
     el(s,'rect',{x:kx,y:ky+i*11,width:7,height:7,rx:1.5,fill:col});
     txt(s,{x:kx+11,y:ky+i*11+6.5,'font-size':7.5,'font-weight':600,fill:c.muted},zh()?(PROV_ZH[p.label]||p.label):p.label)});}
  const mx=X(me[0]),my=Y(me[1]);
  el(s,'polygon',{points:`${mx},${my-7.8} ${mx+7.8},${my} ${mx},${my+7.8} ${mx-7.8},${my}`,fill:c.hero,stroke:c.bg,'stroke-width':1.5,class:'pop'});
  txt(s,{x:mx+11,y:my+4,'font-size':10,'font-weight':800,fill:c.ink,stroke:c.bg,'stroke-width':3,'paint-order':'stroke'},NAME());
  txt(s,{x:x0+W/2,y:y0+H+22,'text-anchor':'middle','font-size':8,'font-weight':700,fill:c.muted,'letter-spacing':'.1em'},'PC1');
  txt(s,{x:x0-30,y:y0+H/2,'text-anchor':'middle','font-size':8,'font-weight':700,fill:c.muted,'letter-spacing':'.1em',transform:`rotate(-90 ${x0-30} ${y0+H/2})`},'PC2');
  foot(s,c,520,462,zh()?'空心方块 = 古代个体（投影）· 圆点 = 现代个体 · 彩色圆 = 汉族样本按省份聚合（色块见左上）· 菱形 = '+NAME():'open squares = ancient (projected) · dots = present-day · coloured circles = Han samples grouped by province (key at top left) · diamond = '+NAME());
});

/* ── distance vs time ── */
/* ── AN6: 离线地图的共用坐标变换、年代筛选与绘制 ─────────────────────────────
   底图是世界范围的等距圆柱投影、固定 viewBox 0 0 360 180（每度 1 单位，x = lon+180,
   y = 90-lat）。地图上每个点都必须走同一个变换，否则点会与轮廓错位。不引入 GIS 库或
   在线瓦片：底图是一次性转好的静态 path，内嵌在页面里。 */

const GEO_VB = { w: 360, h: 180 };
const geoXY = (lat, lon) => ({ x: (Number(lon) + 180) * (GEO_VB.w / 360),
                               y: (90 - Number(lat)) * (GEO_VB.h / 180) });
/* 注意 Number(null) === 0、Number('') === 0：直接 Number() 会把"没有坐标"变成几内亚湾的 (0,0)。
   §7 明确要求空值不能变成 (0,0)，所以先排除 null/undefined/空串，再判断数值与范围。 */
const geoValid = (lat, lon) => {
  if (lat === null || lat === undefined || lat === '' || lon === null || lon === undefined || lon === '') return false;
  const a = Number(lat), b = Number(lon);
  return Number.isFinite(a) && Number.isFinite(b) && Math.abs(a) <= 90 && Math.abs(b) <= 180;
};

/* §7 的筛选约定：判断记录的时间区间与所选范围是否相交。只有均值者按点处理；区间未知的
   记录不猜——它既不进筛选结果，也不被当成"远古"或"现代"。 */
function overlapsAge(minBP, maxBP, loBP, hiBP) {
  return Number.isFinite(minBP) && Number.isFinite(maxBP) && maxBP >= loBP && minBP <= hiBP;
}
function ageInRange(rec, loBP, hiBP) {
  const lo = Number(rec && rec.date_min_bp), hi = Number(rec && rec.date_max_bp);
  const mean = Number(rec && rec.date_mean_bp);
  if (Number.isFinite(lo) && Number.isFinite(hi)) return overlapsAge(lo, hi, loBP, hiBP);
  if (Number.isFinite(mean)) return overlapsAge(mean, mean, loBP, hiBP);
  return null;                       // 年代未知
}

/* AN6 复审：地图行的选取提成纯函数，卡片与回归测试共用同一份（不再各写一套）。
   - 只画合格记录：eligible===false 是 09b/04b 的唯一门槛判定，地图不自行放宽。
   - "全部"范围也显示年代未知的记录：现代参考个体普遍无年代，不能被古代时间轴藏起来；
     具体年代范围内，未知的仍既不算命中也不算排除。 */
function geomapRows(records, loBP, hiBP) {
  const okRows = (records || []).filter(r => r.eligible !== false);
  const inRange = [], outRange = [], noDate = [];
  const isAll = (loBP <= 0 && hiBP >= 1000000);
  okRows.forEach(r => {
    const hit = ageInRange(r, loBP, hiBP);
    if (hit === null) (isAll ? inRange : noDate).push(r);
    else if (hit) inRange.push(r); else outRange.push(r);
  });
  return { inRange, outRange, noDate };
}

/* AN6 复审：09b 的 groups 是 {modern:[...], ancient:[...]}，1000G（04b）的是平铺数组。旧代码
   字典形状只取 .ancient——现代主图（默认视图）按 kind='modern' 过滤后一个不剩，排名/前五名/
   列表全空。归一化把两层的**已算好**结果拼起来：不重算任何排名，各行保留自己梯队里算出的
   rank；并列（两个梯队的第 1 名）按 distance_mean 稳定排序，两边的第 1 名都可见。 */
function analysisGroups(a){
  const g=(a||{}).groups;
  if(Array.isArray(g)) return g;
  if(g&&typeof g==='object') return [].concat(g.modern||[],g.ancient||[]);
  return [];
}

/* 视图分层（AN6 批一，纯函数与测试同源）：modern/ancient 只看各自 kind 的记录，'all' 保持
   原样（含未分层记录）。现代群体是默认主视图——现代参考个体没有年代，古代时间轴属于古代视图，
   两个视图不该被同一个"全部"混在一起。 */
function kindView(records, kind) {
  if (kind !== 'modern' && kind !== 'ancient') return records || [];
  return (records || []).filter(r => String(r.kind || '') === kind);
}

/* 地图取景（纯函数，与回归测试同源）：点全部集中在某一区域（本样本全在东亚）时，整幅世界
   底图把内容压到角落——视口改为跟随当前已落点的包围盒，按其尺寸留边距；跨度下限（70×45
   经纬度）保住区域上下文，不至于一两个点就放大到街区级。无有效坐标、或点横跨全球时回退
   全世界。不做经度折叠（本数据集没有跨 180° 的样本），视口始终夹在世界边界内。 */
function geomapViewport(rows, minW, minH, padFrac) {
  minW = minW || 70; minH = minH || 45; padFrac = padFrac == null ? 0.25 : padFrac;
  const pts = (rows || []).filter(r => geoValid(r.latitude, r.longitude))
                          .map(r => geoXY(r.latitude, r.longitude));
  if (!pts.length) return { x0: 0, y0: 0, w: 360, h: 180 };
  const x0 = Math.min(...pts.map(p => p.x)), x1 = Math.max(...pts.map(p => p.x));
  const y0 = Math.min(...pts.map(p => p.y)), y1 = Math.max(...pts.map(p => p.y));
  const mx = Math.max(4, padFrac * (x1 - x0)), my = Math.max(4, padFrac * (y1 - y0));
  const w = Math.min(360, Math.max(minW, (x1 - x0) + 2 * mx));
  const h = Math.min(180, Math.max(minH, (y1 - y0) + 2 * my));
  const cx = (x0 + x1) / 2, cy = (y0 + y1) / 2;
  return { x0: Math.max(0, Math.min(360 - w, cx - w / 2)),
           y0: Math.max(0, Math.min(180 - h, cy - h / 2)), w, h };
}

/* 共用地图：底图 + 点 + 选择联动。现代图与古代图只是传入的 locations 不同——同一段绘制逻辑，
   避免两个视图对"精度""无坐标"给出不一致的处理。返回 { placed, unplaced } 计数。 */
/* AN6 复审：drawGeoMap 自引入以来没有任何调用方（geomap 卡片自绘），却留着一套与卡片不同的
   精度/距离处理——"辅助函数留着、实际渲染另写一套"正是复审点名的问题。删除；卡片即唯一实现。 */

/* ── AN6 地理分布：把合格记录按采样/发现地点画在世界底图上，并可按年代范围查看 ─────────
   位置是**参考样本的来源地**，不是把目标样本投成某个坐标（§7）。筛选只是换一个查看已算好的
   结果：不重建 PCA，也不在前端重算任何排名。时间轴左古右今，同时给出 BP 与公元（BP 基准 1950）。*/
reveal('geomap',(s,c)=>{
  // 复审 AN6-P1：原先固定优先 AADR，不看数据说哪一份是默认。现在按 default_analysis_id 选；它为空
  // （没有任何分析处于 ok）时不画，而不是拿另一份数据集的位置冒充当前结果。
  const A=(D.ancestry&&D.ancestry.analyses)||[], Z=zh();
  const _defId=(D.ancestry&&D.ancestry.default_analysis_id)||'';
  const a=(_defId&&A.find(x=>String(x.analysis_id)===String(_defId)))|| (_defId?null:A[0]) || {};
  const groups=analysisGroups(a);
  const rankedAll=groups.filter(g=>g.rank&&!g.small_group)
    .sort((x,y)=>(x.rank-y.rank)||(Number(x.distance_mean||0)-Number(y.distance_mean||0)));
  const VB={w:900,h:470}, pad={l:26,t:26}, mapH=300;
  const plotW=VB.w-2*pad.l, plotH=mapH-2*pad.t;
  // sc/ox/oy 随 draw() 里的取景变（见下）；初值=全世界，内容为空时兜底。
  let sc=Math.min(plotW/360,plotH/180);
  let ox=pad.l+(plotW-360*sc)/2, oy=pad.t+(plotH-180*sc)/2;
  const P=r=>{const q=geoXY(r.latitude,r.longitude);return {x:ox+q.x*sc, y:oy+q.y*sc};};
  // 年代范围（BP）。0 = 现在，越往左越古老；1950 是 BP 基准。
  const RANGES=[{l:0,h:1000000,zh:'全部',en:'all'},
                {l:0,h:1500,zh:'1500 BP 以内',en:'< 1500 BP'},
                {l:1500,h:5000,zh:'1500–5000 BP',en:'1500–5000 BP'},
                {l:5000,h:1000000,zh:'5000 BP 以上',en:'> 5000 BP'}];
  // 视图分层（AN6 批一）：现代群体是默认主视图；古代视图自带时间轴；'all' 是原来的混合视图。
  // 切换只是换一批已算好的记录/分组查看，不重算任何统计。
  const KINDS=[{k:'modern',zh:'现代群体',en:'modern'},{k:'ancient',zh:'古代',en:'ancient'},{k:'all',zh:'全部',en:'all'}];
  let kv='modern', cur=0, sel=null;    // sel = 选中的记录 ID（稳定键），点与列表共用
  const LBOX=document.getElementById('geomap_list');
  const bandY=mapH+40;
  const draw=()=>{
    clearEl(s);
    const R=RANGES[cur];
    // 选取走 geomapRows/kindView（与测试同一份纯函数）：不合格记录不进视图；"全部"下无年代者
    // 也显示——现代参考个体普遍没有年代，不能被古代时间轴藏起来。
    const { inRange, outRange, noDate } = geomapRows(kindView(a.records, kv), R.l, R.h);
    // 分层视图下排名/前五名/列表也只看该层：现代主图标现代群体，古代图标古代群体。
    const ranked=(kv==='all'?rankedAll:rankedAll.filter(g=>String(g.kind||'')===kv));
    const top5=ranked.slice(0,5).map(g=>String(g.label));
    const unplaced=inRange.filter(r=>!geoValid(r.latitude,r.longitude)).length;
    const placed=inRange.filter(r=>geoValid(r.latitude,r.longitude));
    // 取景随内容（geomapViewport，与测试同源）：视口装下当前全部已落点并留边距，点因此必在
    // 绘图区内；放大后底图超出绘图区的部分 clip 掉。视口即该视图内容的诚实范围——不画没数据
    // 的大洲，也没有把点投出地图。
    const vp=geomapViewport(placed);
    sc=Math.min(plotW/vp.w,plotH/vp.h);
    ox=pad.l+(plotW-vp.w*sc)/2-vp.x0*sc;
    oy=pad.t+(plotH-vp.h*sc)/2-vp.y0*sc;
    const _defs=el(s,'defs',{}), _cp=el(_defs,'clipPath',{id:'geomap_plot'});
    el(_cp,'rect',{x:pad.l,y:pad.t,width:plotW,height:plotH});
    // 取景放大后 ox/oy 常为负（把世界左移出画布）。use 的 x/y 定位+clip 在个别引擎里不被应用，
    // 底图整体落回原点而点仍在放大位置（人工复核第 3 条截图：北美居左、点悬大西洋）。定位与裁剪
    // 改挂到 g 上、use 只用 transform 平移缩放——transform 是所有引擎一致执行的定位方式，几何等价。
    const _mg=el(s,'g',{'clip-path':'url(#geomap_plot)'});
    el(_mg,'use',{href:'#world_land',x:0,y:0,width:360,height:180,
      transform:`translate(${ox} ${oy}) scale(${sc})`,fill:c.grid,'fill-opacity':.55,stroke:'none'});
    const dists=placed.map(r=>Number(r.distance_to_target)).filter(Number.isFinite);
    const dmax=dists.length?Math.max(...dists):0;
    placed.forEach(r=>{
      const p=P(r), lab=String(r.source_population_id||r.label||''), rid=String(r.record_id||'');
      const isSel=sel&&String(r.record_id||'')===sel;
      const isTop=top5.includes(lab)||isSel;
      // 精度未知不冒充遗址级（AN6）：site/region 是元数据给的判定，缺了就是未知——形状仍可画，
      // 但提示与图例如实说"未知"，不因"有坐标"就当成 site。
      const prec=r.location_precision, site=prec==='site', region=prec==='region';
      const precTxt=site?(Z?'遗址级':'site'):(region?(Z?'地区级':'region'):(Z?'精度未知':'precision unknown'));
      const t=dmax>0&&Number.isFinite(Number(r.distance_to_target))?Math.min(1,Number(r.distance_to_target)/dmax):1;
      const op=isSel?1:(isTop?1:(0.18+0.5*(1-t)));
      const n=region?el(s,'rect',{x:p.x-(isTop?3.4:2),y:p.y-(isTop?3.4:2),width:isTop?6.8:4,height:isTop?6.8:4,rx:1,fill:isTop?c.hero:c.data,'fill-opacity':op,stroke:c.bg,'stroke-width':.4,class:'pop'})
                  :el(s,'circle',{cx:p.x,cy:p.y,r:isTop?4:2.1,fill:isTop?c.hero:c.data,'fill-opacity':op,stroke:c.bg,'stroke-width':.4,class:'pop'});
      n.setAttribute('tabindex','0'); n.setAttribute('role','button');
      const pick=()=>{ sel=(sel===rid)?null:rid; draw(); };
      n.addEventListener('click',pick);
      n.addEventListener('keydown',(e)=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();pick()}});
      tip(n,`${lab||r.record_id} · ${r.locality||''} ${precTxt}`+
            (Number.isFinite(Number(r.date_mean_bp))?` · ${fmt(Math.round(r.date_mean_bp))} BP`:'')+
            (Number.isFinite(Number(r.distance_to_target))?` · d=${Number(r.distance_to_target).toFixed(4)}`:''));
    });
    // 前五名标注：按 member_ids+location_id 找代表点（AN6）——同组多遗址时按 label 找"首个点"
    // 会标到别的遗址上；member_ids 缺席的旧数据才退回 label 匹配。
    ranked.slice(0,5).forEach((g,i)=>{
      const ids=new Set((g.member_ids||[]).map(String));
      const hit=placed.find(r=>ids.has(String(r.record_id))&&String(r.location_id||'')===String(g.location_id||''))
             || placed.find(r=>ids.has(String(r.record_id)))
             || placed.find(r=>String(r.source_population_id)===String(g.label));
      if(!hit) return;
      const p=P(hit);
      txt(s,{x:p.x+6,y:p.y-3-i*9,'font-size':8,'font-weight':700,fill:c.ink,stroke:c.bg,'stroke-width':2.2,'paint-order':'stroke'},
          `${g.rank}. ${ANC_ZH[g.label]||g.label} · n=${g.n}`);
    });
    txt(s,{x:pad.l,y:16,'font-size':8,'font-weight':600,fill:c.faint,'letter-spacing':'.06em'},
        Z?'参考样本的来源地（不是目标样本的坐标）':'WHERE THE REFERENCE SAMPLES COME FROM (NOT THE TARGET)');
    // ── 视图选择（AN6 批一）：现代群体主图 / 古代 / 全部。可聚焦 + Enter；只换查看层，不重算。
    const bx0=pad.l+40, bx1=VB.w-pad.l-40, by=bandY;
    KINDS.forEach((v,i)=>{
      const x0=bx0+i*((bx1-bx0)/KINDS.length), w=(bx1-bx0)/KINDS.length-8, y=by+28;
      const on=v.k===kv;
      const b=el(s,'rect',{x:x0,y:y,width:w,height:18,rx:2,fill:on?c.hero:'none',stroke:on?c.hero:c.grid,'stroke-width':on?0:.8,
                           class:'fade'});
      b.setAttribute('tabindex','0'); b.setAttribute('role','button');
      const go=()=>{ if(kv!==v.k){kv=v.k;draw();} };
      b.addEventListener('click',go);
      b.addEventListener('keydown',(e)=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();go()}});
      txt(s,{x:x0+w/2,y:y+12.5,'text-anchor':'middle','font-size':8,'font-weight':700,fill:on?c.bg:c.ink},Z?v.zh:v.en);
    });
    if(kv!=='modern'){
      // ── 时间带：左古右今；BP 与公元并列（BP 基准 1950）。现代参考个体没有年代——时间轴与
      // 年代范围是古代/全部视图的事，现代主图不带（AN6 批一：古代时间筛选是独立视图）。
      el(s,'line',{x1:bx0,y1:by,x2:bx1,y2:by,stroke:c.grid,'stroke-width':1.2});
      const X=v=>bx1-(bx1-bx0)*Math.min(v,8000)/8000;      // 8000 BP 以上折到左端
      [0,2000,4000,6000,8000].forEach(v=>{
        el(s,'line',{x1:X(v),y1:by-4,x2:X(v),y2:by+4,stroke:c.faint,'stroke-width':.8});
        // BP → 公元：1950 - BP（BP 以 1950 为基准，2000 BP 是公元前 50 年，不是公元 50 年）。
        // 两套基准写错方向就会把"前 2050 年"说成"公元 2050 年"。0 BP 就是 1950 年本身——
        // 标成"今"会把基准年当成当前年份（AN6）。
        const ce=1950-v;
        txt(s,{x:X(v),y:by+15,'text-anchor':'middle','font-size':7.5,fill:c.faint},
            v===0?(Z?'1950（BP 基准）':'1950 (BP base)')
                 :`${fmt(v)} BP`+(ce<=0?`（${Z?'约前':'≈'}${fmt(-ce)}${Z?' 年':''}）`:`（${Z?'约':'≈'}${fmt(ce)}${Z?' 年':' CE'}）`));
      });
      txt(s,{x:bx0-6,y:by+4,'text-anchor':'end','font-size':7.5,'font-weight':700,fill:c.muted},Z?'古老':'older');
      txt(s,{x:bx1+6,y:by+4,'font-size':7.5,'font-weight':700,fill:c.muted},Z?'现代':'recent');
      // 范围选择：可聚焦 + Enter，键盘可用；点击只改查看范围，不重算任何统计
      RANGES.forEach((r,i)=>{
        const x0=bx0+i*((bx1-bx0)/RANGES.length), w=(bx1-bx0)/RANGES.length-8, y=by+50;
        const on=i===cur;
        const b=el(s,'rect',{x:x0,y:y,width:w,height:18,rx:2,fill:on?c.hero:'none',stroke:on?c.hero:c.grid,'stroke-width':on?0:.8,
                             class:'fade'});
        b.setAttribute('tabindex','0'); b.setAttribute('role','button');
        const go=()=>{ if(cur!==i){cur=i;draw();} };
        b.addEventListener('click',go);
        b.addEventListener('keydown',(e)=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();go()}});
        txt(s,{x:x0+w/2,y:y+12.5,'text-anchor':'middle','font-size':8,'font-weight':700,fill:on?c.bg:c.ink},Z?r.zh:r.en);
      });
    }
    // ── 列表：与地图共用同一个 sel（稳定 ID），两边互选；缺坐标的行标出来但仍可选中
    if(LBOX){
      LBOX.textContent='';
      ranked.slice(0,12).forEach(g=>{
        // 与地图共用同一个 sel（稳定 ID），两边互选；行的定位点按 member_ids+location_id 找
        // （AN6）——同组多遗址时按 label 取首个会把选中带到别的遗址。member_ids 缺席才退回 label。
        const ids=new Set((g.member_ids||[]).map(String));
        const hits=placed.filter(r=>ids.has(String(r.record_id))&&String(r.location_id||'')===String(g.location_id||''));
        const fallback=!hits.length?placed.filter(r=>String(r.source_population_id)===String(g.label)):[];
        const cand=hits.length?hits:fallback;
        const rid=cand.length?String(cand[0].record_id||''):'';
        const row=document.createElement('div');
        row.tabIndex=0; row.setAttribute('role','button');
        row.style.cssText='cursor:pointer;padding:1px 3px;border-radius:3px'+(rid&&rid===sel?';background:currentColor;opacity:.14':'');
        row.textContent=`${g.rank}. ${ANC_ZH[g.label]||g.label} · n=${g.n}`+
          (rid?` · d=${Number(g.distance_mean).toFixed(4)}`:(Z?' · 未定位':' · unplaced'));
        const pick=()=>{ if(!rid) return; sel=(sel===rid)?null:rid; draw(); };
        row.addEventListener('click',pick);
        row.addEventListener('keydown',(e)=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();pick()}});
        LBOX.appendChild(row);
      });
    }
    // 计数：分别说明"范围外"与"年代未知"，后者既不算命中也不算排除；"全部"下无年代者照常显示。
    // 现代视图没有年代维度，计数照给但时间轴说明换成视图说明（AN6 批一）。
    const note=document.getElementById('geomap_note');
    const _kv=KINDS.find(v=>v.k===kv)||KINDS[0];
    if(note) note.textContent=(Z
      ? `视图：${_kv.zh}`+(kv==='modern'?'（现代参考无年代，时间轴与年代范围见"古代/全部"视图）':` · 范围：${R.zh}`)
        +` · 命中 ${inRange.length} 条（其中 ${placed.length} 条已定位、${unplaced} 条缺坐标不落点）· 范围外 ${outRange.length} 条 · 年代未知 ${noDate.length} 条（"全部"下照常显示，其余范围不计入筛选）· 深色 = 距离最近的前五名 · 圆点 = 遗址级或精度未知，方框 = 地区级 · 底图 Natural Earth 1:110m（public domain）`
      : `view: ${_kv.en}`+(kv==='modern'?' (modern references carry no date; see the ancient/all views for the timeline)':` · range: ${R.en}`)
        +` · ${inRange.length} records (${placed.length} placed, ${unplaced} without coordinates) · ${outRange.length} outside · ${noDate.length} with an unknown date (shown under "all", excluded from ranged filtering) · emphasised = five closest · circles = site or unknown precision, squares = region-level · base map Natural Earth 1:110m (public domain)`)
      .replace(/\*\*/g,'');
    foot(s,c,VB.w,VB.h-6,Z?'筛选只是换个范围查看已算好的结果，不重建 PCA、也不重算排名 · 底图自动放大到当前视图有内容的区域 · 离线底图，无外部请求'
                          :'filtering only changes which computed records are shown; no PCA rebuild, no re-ranking · the base map auto-zooms to where the current view has content · offline base map, no external requests');
  };
  draw();
});

/* ── affinity ranking (W-T1). Replaces the distance-vs-age scatter: a group mean and a single
      genome's minimum are different statistics, so the old figure mixed them on one y-axis with a
      time axis that carried no signal (r=0.09). Here both kinds share one distance axis, sorted,
      with n and date per row, and the closest genome is shown against its own group below. ── */
reveal('timedist',(s,c)=>{
  const rows=D.ho_affinity||[], strip=D.ho_affinity_strip||[], tg=D.ho_target_group, Z=zh(), a=AC();
  const X0=152, X1=300, top=38, rh=13;
  if(!rows.length){txt(s,{x:12,y:22,'font-size':10,fill:c.faint},Z?'本版没有古 DNA 投影表（09b 未产出）。':'No ancient-DNA projection table in this build (step 09b produced none).');return}
  const nm=r=>Z?(r.name_zh||ANC_ZH[r.label]||r.label.replace(/_/g,' ')):r.label.replace(/_/g,' ');
  const mods=rows.filter(r=>r.kind==='modern'), anc=rows.filter(r=>r.kind==='ancient');
  const ancTop=anc.slice(0,20), keep=new Set(mods.concat(ancTop).map(r=>r.label));
  const shown=rows.filter(r=>keep.has(r.label)), rest=anc.slice(20);   // rows is already sorted by d
  const dmax=Math.max(...shown.map(r=>r.d))*1.10, X=v=>X0+(X1-X0)*v/dmax;
  const COL=r=>r.kind==='modern'?c.data2:(r.region==='north'?a.n:(r.region==='south'?a.s:c.faint));
  const yEnd=top+shown.length*rh;
  [[Z?'现代':'modern',c.data2,0],[Z?'古代·北方':'ancient north',a.n,1],
   [Z?'古代·南方':'ancient south',a.s,2],[Z?'未分类':'unclassified',c.faint,3]].forEach(l=>{
    el(s,'rect',{x:98+l[2]*82,y:8,width:8,height:7,rx:1,fill:l[1]});
    txt(s,{x:110+l[2]*82,y:14,'font-size':7.4,fill:c.muted},l[0])});
  [0,0.02,0.04,0.06].forEach(v=>{if(v>dmax)return;
    el(s,'line',{x1:X(v),y1:top-6,x2:X(v),y2:yEnd,stroke:c.grid,'stroke-width':.5});
    txt(s,{x:X(v),y:yEnd+11,'text-anchor':'middle','font-size':7.4,fill:c.faint},v.toFixed(2))});
  const m1=mods[0];
  if(m1){el(s,'line',{x1:X(m1.d),y1:top-6,x2:X(m1.d),y2:yEnd,stroke:c.ink,'stroke-width':1,'stroke-dasharray':'4 3'});
    txt(s,{x:X(m1.d)+4,y:top-10,'font-size':7.4,'font-weight':700,fill:c.ink},
      (Z?'最近现代人群（':'nearest modern (')+nm(m1)+(Z?'）':'')+' d='+m1.d.toFixed(4));}
  shown.forEach((r,i)=>{const y=top+i*rh, col=COL(r);
    txt(s,{x:X0-8,y:y+3,'text-anchor':'end','font-size':7.6,'font-weight':r.kind==='modern'?700:500,fill:c.ink},nm(r));
    const t=el(s,'line',{x1:X0,y1:y,x2:Math.max(X(r.d),X0+1.5),y2:y,stroke:col,'stroke-width':2.4,class:'fade',style:`animation-delay:${i*.012}s`});
    tip(t,`${r.label} · ${r.kind} · n=${r.n} · d=${r.d.toFixed(5)}`+(r.date_mean?` · ${fmt(r.date_mean)} BP`:''));
    el(s,'circle',{cx:X(r.d),cy:y,r:2.6,fill:col});
    txt(s,{x:X1+6,y:y+3,'font-size':7,fill:c.muted},`n=${r.n}`+(r.date_mean?` · ≈${fmt(r.date_mean)} BP`:''));
  });
  if(rest.length){const y=yEnd+5;
    txt(s,{x:X0-8,y:y+3,'text-anchor':'end','font-size':7.4,fill:c.faint},(Z?'其余 ':'rest of ')+rest.length+(Z?' 个古代群体':' ancient groups'));
    el(s,'line',{x1:X0,y1:y,x2:X(rest[0].d),y2:y,stroke:c.grid,'stroke-width':2.4,'stroke-dasharray':'2 2'});
    txt(s,{x:X1+6,y:y+3,'font-size':7,fill:c.faint},'d ≥ '+rest[0].d.toFixed(3));}
  // The strip: the closest groups' own members on the same axis, so the gap between the closest
  // genome and its group mean is visible rather than asserted in a caption.
  const sy0=yEnd+(rest.length?24:10)+30, groups=strip.concat((tg&&!tg.in_strip)?[tg]:[]);
  txt(s,{x:0,y:sy0-12,'font-size':8,'font-weight':600,fill:c.faint,'letter-spacing':'.06em'},
      Z?'最近的古代群体：组内个体（竖线 = 组均值）':'CLOSEST ANCIENT GROUPS: MEMBERS (TICK = GROUP MEAN)');
  groups.forEach((g,i)=>{const y=sy0+i*24;
    txt(s,{x:X0-8,y:y+3,'text-anchor':'end','font-size':7.8,'font-weight':700,fill:c.ink},nm(g));
    el(s,'line',{x1:X0,y1:y,x2:X1+30,y2:y,stroke:c.grid,'stroke-width':.5});
    el(s,'line',{x1:X(g.mean_d),y1:y-7,x2:X(g.mean_d),y2:y+7,stroke:c.ink,'stroke-width':1.2});
    g.members.forEach(m=>{
      // 高亮的是**数据里那个最近个体**，不是写死的某个 ID：换样本时它自动跟着换。
      const hi=!!((D.ho_near_individual||[])[0] && String(m.iid)===String(D.ho_near_individual[0].iid) && tg && tg.label===g.label);
      const dot=el(s,'circle',{cx:X(m.d),cy:y,r:hi?4.2:3,fill:hi?c.hero:COL(g),stroke:c.bg,'stroke-width':1.1,class:'pop'});
      tip(dot,`${m.iid} · d=${m.d.toFixed(5)} · call_rate ${m.call_rate} · ${fmt(m.date)} BP`);
      if(hi)txt(s,{x:X(m.d),y:y-10,'text-anchor':'middle','font-size':7.4,'font-weight':800,fill:c.hero,stroke:c.bg,'stroke-width':2.4,'paint-order':'stroke'},
        String(m.iid).replace(/\.(SG|AG)$/,'')+' · call_rate '+m.call_rate.toFixed(2));});
    txt(s,{x:X1+36,y:y+3,'font-size':7,fill:c.muted},(Z?'均值 ':'mean ')+g.mean_d.toFixed(4)+' · n='+g.n);});
  const ny=sy0+groups.length*24+14;
  [Z?'距离 = HO-PCA 前 4 主成分的欧氏距离，越小越近；现代与古代在同一把尺上。'
      :'distance = Euclidean distance in the first 4 HO-PCA components, smaller = closer; modern and ancient share one scale',
   Z?`单个体的距离噪声大：${((D.ho_near_individual||[])[0]||{}).iid||'该个体'} 与其群体均值的差距不构成亲缘结论。`
      :`a single genome is noisy: the gap between ${((D.ho_near_individual||[])[0]||{}).iid||'that individual'} and its group mean is not a kinship result`,
   Z?'本图（全局 PCA）与 FLARE 局部祖源（南北面板）是不同尺度，互不替代、也不矛盾。'
      :'this figure (global PCA) and the FLARE local-ancestry panels are different scales: neither replaces nor contradicts the other',
   Z?'区域着色依 panel/aadr_site_regions.tsv（秦岭—淮河分界）；未分类不猜色。'
      :'regions follow panel/aadr_site_regions.tsv (Qinling-Huaihe line); unclassified rows are not colour-guessed'
  ].forEach((line,i)=>txt(s,{x:0,y:ny+i*11,'font-size':7.4,fill:c.faint},line));
  foot(s,c,420,ny+4*11+12,Z?'条形 = 到样本的距离 · 虚线 = 最近现代人群 · 竖线 = 组均值':'bar = distance to the sample · dashed = nearest modern group · tick = group mean');
});
/* ── f3 statistics card (28_f3_stats.py). Built like the lineage card (JS-filled div), not a
      reveal() figure: the copy is bilingual inline and the whole card hides when D.f3 is null,
      which the sections table already explains. Numbers only; no group names are hardcoded --
      sets, labels and n all come from the payload. ── */
/* Pure helpers for the f3 card (node-tested in tests/test_f3_card.py). The corrected
   producer (28) emits null f3/se/z when a statistic is unavailable (zero-SE, non-finite),
   so every displayed number and verdict comes from the row itself -- Number(null)===0 means
   the old _sg printed "+0.0" for unavailable rows, and the old note asserted "all positive"
   about the current sample instead of reading the data. */
function f3Num(v,d){ if(v==null||!isFinite(v)) return null; return (v>0?'+':'')+Number(v).toFixed(d); }
function f3Verdict(r){
  const fin=(v)=>v!=null&&isFinite(v);   // isFinite(null)===true (coerced to 0) -- guard explicitly
  if(r==null||!fin(r.f3)||!fin(r.se)||!fin(r.z)) return 'unavailable';
  return (r.f3<0&&Math.abs(r.z)>3)?'sig-neg':'ns';
}
function f3Cell(r){ const f=f3Num(r.f3,4),s=f3Num(r.se,4),z=f3Num(r.z,1);
  if(f==null||s==null) return '—';
  return f+' ±'+s.replace('+','')+'  (Z='+(z==null?'—':z)+')'; }
function f3Note(adm){
  const av=adm.filter(r=>r.f3!=null&&r.se!=null&&isFinite(r.f3)&&isFinite(r.se));
  if(!av.length) return adm.length?'all-unavailable':'none';
  if(av.some(r=>f3Verdict(r)==='sig-neg')) return 'has-sig-neg';
  return 'no-sig-neg';
}
const renderF3=()=>{
  const box=document.getElementById('f3card'); if(!box) return;
  const FD=D.f3; if(!FD){box.style.display='none';return}
  const M=FD.modern||null, A=FD.ancient||null, Z=zh();
  const rows=(M&&M.outgroup_f3)||[];
  if(!rows.length&&!(A&&(A.admixture||[]).length)){box.style.display='none';return}
  const _sg=(v,d)=>{const s=f3Num(v,d);return s==null?'—':s;};
  // Ranking bars: the spread between groups is a tiny fraction of the absolute value, so the
  // axis is floored near the data (lo = min - 3*SE); the whiskers carry the real uncertainty.
  // Rows whose f3/SE are unavailable (null) draw no bar -- they get an explicit 不可用 marker.
  const x0=150,X1=310,top=30,rh=19;
  const draw=rows.filter(r=>Number.isFinite(r.f3)&&Number.isFinite(r.se));
  const lo=draw.length?Math.min(...draw.map(r=>r.f3-3*r.se)):0, hi=draw.length?Math.max(...draw.map(r=>r.f3+3*r.se)):1;
  const X=v=>x0+(X1-x0)*(v-lo)/(hi-lo);
  const adm=(M&&M.admixture||[]).map(r=>({...r,scope:'m'})).concat((A&&A.admixture||[]).map(r=>({...r,scope:'a'})));
  const gM=(M&&M.groups)||{};
  const H=top+(rows.length*rh||0)+(rows.length?16:0)+(adm.length?adm.length*18+26:0)+16;
  box.innerHTML=
    '<div class="claim" style="margin-top:34px">'+(Z?'f3 统计：与各参照群的共享漂移':'f3 statistics: shared drift with each reference group')+'</div>'
   +'<div class="sub">'+(Z
      ?'外群 f3 越大 = 与该群共享的漂移越多（配对差异的 Z 是排序的正式依据）；混合 f3 显著为负才是两群混合的证据'
      :'larger outgroup-f3 = more shared drift with that group (adjacent-rank contrast Z values are the formal basis for the ordering); a significantly negative admixture-f3 is the signature of two-source admixture')+'</div>'
   +'<div class="fig"><svg id="f3" viewBox="0 0 900 '+Math.max(H,140)+'" preserveAspectRatio="xMidYMid meet"></svg></div>'
   +'<p class="note" id="n_f3" style="max-width:86ch"></p>';
  const s=document.getElementById('f3'), c=C(), a=AC(), note=document.getElementById('n_f3');
  let y=top;
  if(rows.length){
    txt(s,{x:0,y:14,'font-size':8,'font-weight':600,fill:c.faint,'letter-spacing':'.08em'},
      Z?('外群 F3（'+(M.outgroup||'')+'；样本，P） · 共享漂移'):('OUTGROUP-F3('+(M.outgroup||'')+'; TARGET, P) · SHARED DRIFT'));
    [lo,lo+(hi-lo)/2,hi].forEach(v=>{el(s,'line',{x1:X(v),y1:top-6,x2:X(v),y2:top+rows.length*rh-6,stroke:c.grid,'stroke-width':.5});
      txt(s,{x:X(v),y:top+rows.length*rh+8,'text-anchor':'middle','font-size':7.4,fill:c.faint},v.toFixed(4));});
    rows.forEach((r,i)=>{const yy=y+i*rh, col=r.kind==='pool'?c.ink:c.data2;
      txt(s,{x:x0-10,y:yy+3,'text-anchor':'end','font-size':8.4,'font-weight':r.kind==='pool'?700:600,fill:c.ink},
        Z?(r.label_zh||r.set):(r.label_en||r.set));
      if(!Number.isFinite(r.f3)||!Number.isFinite(r.se)){
        txt(s,{x:X1+10,y:yy+3,'font-size':7.6,'font-weight':700,fill:c.muted},
          Z?'不可用（SE=0/非有限）':'unavailable');
        txt(s,{x:X1+62,y:yy+3,'font-size':7,fill:c.muted},'n='+r.n);return}
      el(s,'line',{x1:X(r.f3-r.se),y1:yy,x2:X(r.f3+r.se),y2:yy,stroke:c.faint,'stroke-width':1});
      const b=el(s,'line',{x1:X(Math.max(r.f3-r.se,lo)),y1:yy,x2:X(r.f3),y2:yy,stroke:col,'stroke-width':3,class:'fade',style:`animation-delay:${i*.04}s`});
      tip(b,`${r.set} · f3=${r.f3.toFixed(6)} ±${r.se.toFixed(6)} · Z=${_sg(r.z,1)} · n=${r.n}`);
      el(s,'circle',{cx:X(r.f3),cy:yy,r:2.4,fill:col});
      txt(s,{x:X1+10,y:yy+3,'font-size':7.6,'font-weight':700,fill:c.ink},r.f3.toFixed(5));
      txt(s,{x:X1+62,y:yy+3,'font-size':7,fill:c.muted},'n='+r.n);});
    y+=rows.length*rh+8;
    (M.contrasts||[]).forEach((ct,i)=>{const yy=top+i*rh+rh/2+3;
      txt(s,{x:470,y:yy,'font-size':7.4,'font-weight':Math.abs(ct.z)>=3?800:400,fill:Math.abs(ct.z)>=3?c.ink:c.muted},
        ct.a+' − '+ct.b+':  Δ='+_sg(ct.diff,4)+'  (Z='+_sg(ct.z,1)+')');});
    y+=18;
  }
  if(adm.length){
    txt(s,{x:0,y:y+4,'font-size':8,'font-weight':600,fill:c.faint,'letter-spacing':'.08em'},
      Z?'混合检验 · ADMIXTURE-F3（样本；来源A，来源B）':'ADMIXTURE-F3(TARGET; SOURCE A, SOURCE B)');
    y+=14;
    adm.forEach((r,i)=>{const yy=y+i*18;
      const lab=r.scope==='a'
        ?(Z?(r.label_zh||((r.a+'×'+r.b))):(r.label_en||(r.a+' x '+r.b)))
        :(Z?(((gM[r.a]||{}).label_zh||r.a)+'×'+((gM[r.b]||{}).label_zh||r.b)):(r.a+' x '+r.b));
      txt(s,{x:0,y:yy+3,'font-size':8.4,'font-weight':600,fill:c.ink},
        (r.scope==='a'?(Z?'古参照 · ':'ancient · '):(Z?'现代 · ':'modern · '))+lab);
      txt(s,{x:470,y:yy+3,'font-size':8,'font-weight':700,fill:c.ink},
        f3Cell(r));
      const n=(r.scope==='a')?('n='+r.n_a+'/'+r.n_b):('n='+(((gM[r.a]||{}).n)||'—')+'/'+(((gM[r.b]||{}).n)||'—'));
      txt(s,{x:680,y:yy+3,'font-size':7,fill:c.muted},n);
      const vd=f3Verdict(r);
      tip(el(s,'circle',{cx:462,cy:yy,r:2.4,fill:vd==='sig-neg'?a.s:c.faint}),
        vd==='sig-neg'?(Z?'显著为负：两群混合的证据':'significantly negative: evidence of two-source admixture')
        :vd==='unavailable'?(Z?'不可用：SE=0 或非有限值，无法给出 Z 与结论':'unavailable: zero SE or non-finite value; no Z, no verdict')
                           :(Z?'未检出显著为负信号：不构成两群混合证据（功效受参照分化程度限制）':'no significantly negative signal: not evidence of two-source admixture (power depends on source divergence)'));});
    y+=adm.length*18+6;
  }
  foot(s,c,900,Math.max(H-6,y),
    Z?'条形 = 外群 f3（越右越近） · 须 = ±1 标准误 · 粗体对比 = |Z|≥3 的相邻排名差异'
      :'bar = outgroup-f3 (right = closer) · whisker = ±1 SE · bold contrast = adjacent-rank difference with |Z|≥3');
  const bits=[];
  if(M)bits.push(Z?`现代：${fmt(M.sites)} 个 LD 修剪位点 · ${M.blocks} 个 ${(FD.block_mb||5)}Mb 块 jackknife（1000G）`
                  :`modern: ${fmt(M.sites)} LD-pruned sites · ${M.blocks} ${(FD.block_mb||5)}Mb-block jackknife (1000G)`);
  if(A)bits.push(Z?`古参照：${fmt(A.sites)} 个位点 · ${A.blocks} 块（AADR 区域池，${(A.admixture||[]).map(r=>r.label_zh+' n='+r.n_a+'/'+r.n_b).join('，')}）`
                  :`ancient: ${fmt(A.sites)} sites · ${A.blocks} blocks (AADR region pools, ${(A.admixture||[]).map(r=>r.label_en+' n='+r.n_a+'/'+r.n_b).join(', ')})`);
  // The conclusion sentence is derived from the rows (f3Note), never hardcoded to the current
  // sample's all-positive result; the estimator line discloses what 28 actually computed.
  const verdict=f3Note(adm);
  const concl = verdict==='has-sig-neg'
    ?(Z?'存在显著为负的混合 f3（|Z|>3）：对该对来源构成两群混合信号；幅度解释受参照选择与覆盖限制。'
       :'At least one admixture-f3 is significantly negative (|Z|>3): a two-source admixture signal for that source pair; magnitude interpretation is limited by reference choice and coverage.')
    :verdict==='no-sig-neg'
    ?(Z?'未检出显著为负的混合 f3：这只是"未检出"，不是反证——东亚参照群彼此分化浅，该检验对南北混合的功效本就有限。'
       :'No admixture-f3 is significantly negative: non-detection, not disproof -- East Asian reference groups are weakly differentiated, and the test has little power for north-south admixture.')
    :verdict==='all-unavailable'
    ?(Z?'混合 f3 的标准误为 0 或非有限，无法给出 Z，也就无法给出混合结论。'
       :'Admixture-f3 standard errors are zero or non-finite: no Z, therefore no admixture verdict.')
    :'';
  note.textContent=(bits.join('  ·  ')+(bits.length?'. ':''))+concl
   +(FD.estimator?((concl?'  ':'')+(Z?'估计量：':'estimator: ')+FD.estimator+'  ·  '):'')
   +(Z?'f3 刻画的是等位基因频率的相对接近程度，不据此做族群或籍贯推断。'
      :'f3 measures relative allele-frequency closeness only; no ethnic or geographic-origin inference is drawn from it.');
};

/* ── archaic ideogram ── */
reveal('archaic',(s,c)=>{
  const x0=30,x1=880,rows=11,rh=28,maxL=D.chrlen['1'];
  const by={}; D.archaic.forEach(g=>{(by[g.chrom]=by[g.chrom]||[]).push(g)});
  CHR22.forEach((ch,i)=>{
    const col=i<rows?0:1, row=i%rows, bx=col===0?x0:x0+440, bw=(x1-x0-60)/2;
    const y=20+row*rh, w=bw*D.chrlen[ch]/maxL;
    txt(s,{x:bx-6,y:y+7,'text-anchor':'end','font-size':10,'font-weight':700,fill:c.ink},ch);
    el(s,'rect',{x:bx,y:y,width:w,height:10,rx:2,fill:c.track});
    const cen=D.cen[ch]; if(cen) el(s,'line',{x1:bx+w*cen*1e6/D.chrlen[ch],y1:y-2,x2:bx+w*cen*1e6/D.chrlen[ch],y2:y+12,stroke:c.faint,'stroke-width':1});
    (by[ch]||[]).forEach(g=>{
      const gx=bx+w*g.start/D.chrlen[ch], gw=Math.max(1.5,w*(g.end-g.start)/D.chrlen[ch]);
      const col2=g.src==='Denisovan'?c.data2:AC().a;
      const r=el(s,'rect',{x:gx,y:g.hom?y-1.5:y,width:gw,height:g.hom?13:10,rx:1.5,fill:col2,opacity:g.hom?1:.72});
      tip(r,`chr${ch}:${fmt(g.start)}-${fmt(g.end)} · ${g.src} · ${g.mb} Mb${g.hom?' · homozygous':''}`)});
    txt(s,{x:bx+w+6,y:y+8,'font-size':9,fill:c.muted},(by[ch]||[]).length);
  });
  foot(s,c,900,336,zh()?'绿色 = 尼安德特 · 浅蓝 = 丹尼索瓦 · 加高 = 两条染色体都带 · 右侧 = 片段数':'green = Neanderthal · light blue = Denisovan · taller = on both chromosomes · right = segment count');
});

/* ── mutation spectrum ── */
reveal('spectrum',(s,c)=>{
  const a=AC(), SUB=['C>A','C>G','C>T','T>A','T>C','T>G'];
  const COLS=[c.data2,c.ink,a.s,c.faintdata,a.n,c.data];
  const x0=44,x1=880,y0=30,H=190, data=D.spectrum, max=Math.max(...data.map(d=>d.frac)), bw=(x1-x0)/96, gap=bw*0.22;
  txt(s,{x:x0,y:14,'font-size':8,'font-weight':600,fill:c.faint,'letter-spacing':'.08em'},zh()?'占全部单碱基替换的比例':'SHARE OF ALL SINGLE-BASE SUBSTITUTIONS');
  data.forEach((d,i)=>{const gi=SUB.indexOf(d.sub), x=x0+i*bw, h=H*d.frac/max;
    const r=el(s,'rect',{x:x+gap/2,y:y0+H-h,width:bw-gap,height:h,rx:1,fill:COLS[gi]});
    tip(r,`${d.ctx}: ${(d.frac*100).toFixed(2)}% (${fmt(d.n)})`)});
  el(s,'line',{x1:x0,y1:y0+H,x2:x1,y2:y0+H,stroke:c.faint,'stroke-width':.9});
  SUB.forEach((sb,gi)=>{const xa=x0+gi*16*bw, xb=xa+16*bw;
    el(s,'rect',{x:xa+1,y:y0-12,width:16*bw-2,height:7,rx:1.5,fill:COLS[gi]});
    txt(s,{x:(xa+xb)/2,y:y0+H+16,'text-anchor':'middle','font-size':11,'font-weight':800,fill:c.ink},sb);
    const tot=data.filter(d=>d.sub===sb).reduce((p,q)=>p+q.frac,0);
    txt(s,{x:(xa+xb)/2,y:y0+H+30,'text-anchor':'middle','font-size':9,fill:c.muted},(tot*100).toFixed(1)+'%')});
  const ci=data.map((d,i)=>/\[C>T\]G/.test(d.ctx)?i:-1).filter(i=>i>=0);
  ci.forEach(i=>el(s,'circle',{cx:x0+i*bw+bw/2,cy:y0+H+42,r:2.2,fill:a.s}));
  txt(s,{x:x0+ci[0]*bw-6,y:y0+H+46,'text-anchor':'end','font-size':9,fill:c.muted},zh()?`CpG 位点的 C>T：${F.cpg}%`:`C>T at CpG: ${F.cpg}%`);
  foot(s,c,900,296,zh()?'96 根柱 = 6 种替换 × 前后各一个碱基的 16 种组合':'96 bars = 6 substitution types × 16 flanking-base combinations');
});

/* ── KIR ── */
/* Pure helpers for the KIR card (node-tested in tests/test_kir_card.py). The haplotype call
   and this run's HLA ligand groups come from 24's kir_summary via 30 (D.kir_summary); the
   card used to hardcode the example sample's haplotype and HLA ligand sentence, which
   contradicted this sample's real call. Genes are "unknown" when T1K never ran (empty
   kir.tsv), and only "absent" when the typing actually reports the gene missing. */
function kirRan(rows){ return Array.isArray(rows) && rows.length>0; }
function kirHapLabel(KS,Z){
  const hap=KS&&KS.haplotype;
  if(!hap) return Z?'单倍型：未知（KIR 分型未运行）':'Haplotype: unknown (KIR typing not run)';
  if(/^unavailable/i.test(hap)) return Z?'单倍型：未知（T1K 未运行）':'Haplotype: unknown (T1K not run)';
  if(hap==='AA (two A haplotypes)') return Z?'单倍型：AA（两条 A 单倍型）':'Haplotype: AA (two A haplotypes)';
  if(hap==='Bx (at least one B haplotype)') return Z?'单倍型：Bx（至少一条 B 单倍型）':'Haplotype: Bx (at least one B haplotype)';
  return (Z?'单倍型：':'Haplotype: ')+hap;
}
function kirLigandLine(KS,Z){
  const cg=(KS&&KS.C_ligands)||'', bw=(KS&&KS.B_epitopes)||'', a3=(KS&&KS.A311_ligands)||'';
  const seg=[
    (Z?'HLA-C 配体组：':'HLA-C ligand groups: ')+(cg||'—'),
    (Z?'HLA-B 表位：':'HLA-B epitopes: ')+(bw||'—'),
    (Z?'HLA-A3/A11（KIR3DL2 配体）：':'HLA-A3/A11 (KIR3DL2 ligands): ')+(a3||'—')];
  return seg.join(Z?' · ':'; ');
}
reveal('kir',(s,c)=>{
  const a=AC();
  const ORDER=['KIR3DL3','KIR2DS2','KIR2DL2','KIR2DL3','KIR2DP1','KIR2DL1','KIR3DP1','KIR2DL4','KIR3DL1','KIR3DS1','KIR2DL5A','KIR2DL5B','KIR2DS3','KIR2DS5','KIR2DS1','KIR2DS4','KIR3DL2'];
  const ACT=new Set(['KIR2DS1','KIR2DS2','KIR2DS3','KIR2DS4','KIR2DS5','KIR3DS1']);
  const ran=kirRan(D.kir);
  const map=Object.fromEntries(D.kir.map(k=>[k.gene,k]));
  const x0=30,x1=880,y=64,bw=(x1-x0)/ORDER.length;
  txt(s,{x:x0,y:18,'font-size':8,'font-weight':600,fill:c.faint,'letter-spacing':'.08em'},zh()?'KIR 基因区（19q13.4）按染色体顺序':'THE KIR LOCUS (19q13.4) IN CHROMOSOMAL ORDER');
  el(s,'line',{x1:x0,y1:y+22,x2:x1,y2:y+22,stroke:c.quiet,'stroke-width':1});
  ORDER.forEach((g,i)=>{const k=map[g], on=ran&&k&&k.present, x=x0+i*bw;
    // no typing at all -> unknown (light outline, no fill); typing present -> real absent is dashed
    const r=el(s,'rect',{x:x+2,y:y,width:bw-4,height:44,rx:2.5,fill:on?(ACT.has(g)?a.s:a.n):'none',opacity:on?.9:1,
      stroke:on?'none':(ran?c.faint:c.quiet),'stroke-width':1.2,'stroke-dasharray':on?'':'3 3'});
    tip(r,on?`${g}: ${k.call}`:(ran?`${g}: absent`:`${g}: unknown (KIR typing not run)`));
    txt(s,{x:x+bw/2,y:y-8,'text-anchor':'middle','font-size':8.5,'font-weight':on?700:500,fill:on?c.ink:c.faint,transform:`rotate(-42 ${x+bw/2} ${y-8})`},g.replace('KIR',''));
    if(on) txt(s,{x:x+bw/2,y:y+27,'text-anchor':'middle','font-size':9,fill:c.bg,'font-weight':800},'✓')});
  [[a.n,zh()?'抑制性受体':'inhibitory'],[a.s,zh()?'活化性受体':'activating'],['none',ran?(zh()?'缺失':'absent'):(zh()?'未知（未运行）':'unknown (not run)')]].forEach(([col,lab],i)=>{
    const lx=x0+i*160;
    if(col==='none') el(s,'rect',{x:lx,y:y+44,width:11,height:11,rx:2,fill:'none',stroke:c.faint,'stroke-width':1.2,'stroke-dasharray':'3 3'});
    else el(s,'rect',{x:lx,y:y+44,width:11,height:11,rx:2,fill:col});
    txt(s,{x:lx+16,y:y+53,'font-size':9.5,fill:c.lab},lab)});
  txt(s,{x:x0,y:y+82,'font-size':12,'font-weight':800,fill:c.ink},kirHapLabel(D.kir_summary,zh()));
  txt(s,{x:x0,y:y+99,'font-size':10,fill:c.lab},kirLigandLine(D.kir_summary,zh()));
});

/* ── stat grids and tables ── */
function statGrid(id,rows){const n=document.getElementById(id); if(!n)return;
  n.innerHTML=rows.map(([v,r,l])=>`<div class="kpi"><div class="v">${v}</div><div class="r">${r}</div><div class="l">${l}</div></div>`).join('')}
function deepRender(){
  const hd=document.getElementById('tb_hlad');
  if(hd) hd.innerHTML=D.hla_disease.map(r=>{const st=r.carried, yes=st==='yes', unk=st==='unknown'||st==='';
    const cond=zh()?r.condition:r.condition_en, req=zh()?r.requirement:r.requirement_en, nt=zh()?r.note:r.note_en;
    const cell=unk?(zh()?'未分型':'unknown'):(yes?(zh()?'携带':'yes'):(zh()?'未携带':'no'));
    const col=unk?'#c9a227':(yes?'var(--ancs)':'var(--muted)');
    return `<tr><td class="gene">${cond}</td><td class="rs">${req}</td><td class="g" style="color:${col}">${cell}</td><td>${nt}</td></tr>`}).join('');
  const cd=document.getElementById('tb_cand');
  if(cd) cd.innerHTML=D.candidate.map(r=>`<tr><td class="gene">${gI(r.gene)}</td><td class="rs">${r.rsid} · ${zh()?r.variant:(r.variant_en||r.variant)}</td><td class="g">${r.genotype}</td><td class="rs">${zh()?r.popular_claim:(r.popular_claim_en||'')}</td><td>${zh()?r.what_evidence_supports:(r.evidence_en||'')}</td></tr>`).join('');
  statGrid('k_la',[[`${F.north}<small>%</small>`,zh()?`「${F.la_a}」成分，比第一个校准组平均高 ${F.sdchb} 个标准差`:`share of "${F.la_a}"; ${F.sdchb} sd above the first calibration groupean`,'northern east asian'],
                   [`${F.south}<small>%</small>`,zh()?`「${F.la_b}」成分（校准组 ${100-F.chb}% / ${100-F.chs}%）`:`share of "${F.la_b}" (calibration groups ${100-F.chb}% / ${100-F.chs}% Han ${100-F.chs}%)`,'southern east asian'],
                   [`${F.noise}<small>%</small>`,zh()?'欧洲 + 南亚，方法的噪声底':'European + South Asian, the noise floor','noise floor'],
                   [`${D.la_segments.filter(g=>_P1&&g.anc===_P1).length}`,zh()?`≥0.5 Mb 的「${F.la_b}」片段，最长 ${F.segmax} Mb`:`segments ≥0.5 Mb of "${F.la_b}", longest ${F.segmax} Mb`,'segments']]);
  statGrid('k_arch',[[`${F.archmb}<small>Mb</small>`,zh()?`古老人类片段总长，占常染色体 ${F.archpct}%`:`archaic span, ${F.archpct}% of the autosomes`,'archaic span'],
                     [`${F.archn}`,zh()?`片段数，其中 ${F.archhom} 段为纯合`:`segments, ${F.archhom} of them homozygous`,'segments'],
                     [`${F.neamb}<small>Mb</small>`,zh()?'尼安德特来源':'from Neanderthals','neanderthal'],
                     [`${F.denmb}<small>Mb</small>`,zh()?'丹尼索瓦来源':'from Denisovans','denisovan']]);
  statGrid('k_phase',[[`${F.phpct}<small>%</small>`,zh()?`杂合位点已定相（${F.phased} / ${F.phhet}）`:`heterozygous sites phased (${F.phased} of ${F.phhet})`,'phased'],
                      [`${F.n50}<small>Mb</small>`,zh()?`相位块 N50，共 ${F.blocks} 个块`:`phase block N50, ${F.blocks} blocks`,'block n50'],
                      [`${F.phmax}<small>Mb</small>`,zh()?'最大相位块':'largest phase block','max block']]);
  statGrid('k_som',[[`${F.mtcn}`,zh()?'线粒体拷贝/细胞（正常 100–500）':'mtDNA copies per cell (100–500 normal)','mtdna copy number'],
                    [`${F.ydr}`,zh()?'Y 染色体深度比（<0.45 才提示镶嵌性丢失）':'Y depth ratio (mosaic loss below 0.45)','y depth ratio'],
                    [`${F.chipalt}`,zh()?`${((D.chip_hotspots||{}).n_tested)||0} 个克隆性造血热点共 ${F.chipalt} 条 ALT 读段${(((D.chip_hotspots||{}).hotspots)||[]).length?`（最多：${((D.chip_hotspots.hotspots[0]||{}).hotspot)||''} ${(100*((D.chip_hotspots.hotspots[0]||{}).vaf||0)).toFixed(1)}%）`:''}；低水平读段是噪声敏感的观察，不构成 CHIP 证实，也不排除更低水平克隆`:`${F.chipalt} alt reads across ${((D.chip_hotspots||{}).n_tested)||0} CHIP hotspots${(((D.chip_hotspots||{}).hotspots)||[]).length?` (max: ${((D.chip_hotspots.hotspots[0]||{}).hotspot)||''} ${(100*((D.chip_hotspots.hotspots[0]||{}).vaf||0)).toFixed(1)}%)`:''}; low-level reads are a noise-prone observation, neither confirming CHIP nor excluding lower-level clones`,'clonal haematopoiesis'],
                    [`${F.telk7}`,zh()?`端粒重复读段（共 ${F.teltot} 条读段中），不足以估计端粒长度（k7≥100000 为工作门槛而非验证阈值；绝对长度另需 GC/读长校正）`:`telomeric reads of ${F.teltot} total, too few for a length estimate (k7>=100000 is a working gate, not a validated threshold; absolute lengths would need GC/read-length correction)`,'telomeric reads']]);
  // §7 要求默认视图不重复堆叠个体榜；W-T1 要求 M13 与其群体均值的落差在图（strip）里可见。
  // 两者并存的方式：strip 默认可见，这张逐个体的明细表默认折叠，展开后按同一距离排序。
  const _sum=document.getElementById('ancd_sum');
  if(_sum) _sum.textContent=(zh()?'展开：逐个体的距离明细（与上方 strip 同一距离、同一排序）'
                                :'Expand: per-individual distances (same distances, same order as the strip above)')
    +` · ${(D.ho_near_individual||[]).length}`;
  const tb=document.getElementById('tb_anc');
  // Individual level, not group level: the caption says the nearest ancient individual is M13 of
  // Baiyangcun, and that genome sits inside China_MLBA, whose *group average* ranks far lower. The
  // table used to list group averages only, so the quoted individual appeared nowhere in the report.
  if(tb) tb.innerHTML=(D.ho_near_individual||[]).slice(0,10).map(r=>{const nm=zh()?(ANC_ZH[r.iid]||r.iid):r.iid;
    return `<tr><td class="gene">${nm}</td><td class="rs">${zh()?(ANC_ZH[r.group]||r.group):r.group.replace(/_/g,' ')}</td><td class="num">${fmt(Math.round(r.date))}</td><td class="num">${r.d.toFixed(4)}</td></tr>`}).join('');
  const bb=document.getElementById('tb_blood');
  if(bb) bb.innerHTML=BLOODV3.map(r=>`<tr><td class="gene">${zh()?r.sys_zh:r.sys_en}</td><td class="g">${zh()?r.geno_zh:r.geno_en}</td><td class="gene">${zh()?r.ph_zh:r.ph_en}</td><td>${zh()?r.ev_zh:r.ev_en}</td></tr>`).join('');
}

/* ── circular genome overview ── */
reveal('circos',(s,c)=>{
  const a=AC(), CX=430, CY=430, GAP=0.9*Math.PI/180;
  const total=CHR22.reduce((p,ch)=>p+D.chrlen[ch],0);
  const span=2*Math.PI-CHR22.length*GAP;
  let ang=-Math.PI/2+GAP/2; const arcs={};
  CHR22.forEach(ch=>{const a0=ang,a1=ang+span*D.chrlen[ch]/total; arcs[ch]=[a0,a1]; ang=a1+GAP});
  const P=(r,t)=>[CX+r*Math.cos(t),CY+r*Math.sin(t)];
  const arc=(r0,r1,a0,a1,fill,op)=>{
    const [x0,y0]=P(r1,a0),[x1,y1]=P(r1,a1),[x2,y2]=P(r0,a1),[x3,y3]=P(r0,a0);
    const big=(a1-a0)>Math.PI?1:0;
    return el(s,'path',{d:`M${x0} ${y0}A${r1} ${r1} 0 ${big} 1 ${x1} ${y1}L${x2} ${y2}A${r0} ${r0} 0 ${big} 0 ${x3} ${y3}Z`,fill,opacity:op===undefined?1:op});
  };
  const at=(ch,pos)=>{const A=arcs[ch];if(!A)return 0;const [a0,a1]=A;return a0+(a1-a0)*pos/D.chrlen[ch]};
  const R={lab:410, ideo:[372,392], dens:[318,368], la:[292,314], arch:[272,288], roh:[256,268]};
  CHR22.forEach(ch=>{
    const [a0,a1]=arcs[ch];
    arc(R.ideo[0],R.ideo[1],a0,a1,c.track);
    const cen=D.cen[ch];
    if(cen){const t=at(ch,cen*1e6); const [x0,y0]=P(R.ideo[0],t),[x1,y1]=P(R.ideo[1],t);
      el(s,'line',{x1:x0,y1:y0,x2:x1,y2:y1,stroke:c.faint,'stroke-width':1.4})}
    const am=(a0+a1)/2, [lx,ly]=P(R.lab,am);
    txt(s,{x:lx,y:ly,'text-anchor':'middle','dominant-baseline':'central','font-size':11,'font-weight':700,fill:c.lab},ch);
    for(let p=0;p<D.chrlen[ch];p+=50e6){const t=at(ch,p);const [x0,y0]=P(R.ideo[1],t),[x1,y1]=P(R.ideo[1]+4,t);
      el(s,'line',{x1:x0,y1:y0,x2:x1,y2:y1,stroke:c.faint,'stroke-width':.7})}
  });
  const dens=Object.fromEntries(D.density.map(d=>[d.chrom,d.counts]));
  const allc=D.density.flatMap(d=>d.counts.filter(x=>x>0)).sort((p,q)=>p-q);
  const lo=allc[Math.floor(allc.length*0.02)], hi=allc[Math.floor(allc.length*0.98)];
  el(s,'circle',{cx:CX,cy:CY,r:R.dens[0],fill:'none',stroke:c.grid,'stroke-width':.7});
  el(s,'circle',{cx:CX,cy:CY,r:R.dens[1],fill:'none',stroke:c.grid,'stroke-width':.7});
  CHR22.forEach(ch=>{
    const co=dens[ch]||[]; if(co.length<2) return;
    let d=`M${P(R.dens[0],at(ch,0))[0]} ${P(R.dens[0],at(ch,0))[1]}`;
    co.forEach((v,i)=>{const t=at(ch,Math.min(i*5e6+2.5e6,D.chrlen[ch]));
      const f=Math.max(0,Math.min(1,(v-lo)/(hi-lo)));
      const p=P(R.dens[0]+(R.dens[1]-R.dens[0])*f,t); d+=`L${p[0].toFixed(1)} ${p[1].toFixed(1)}`});
    const [ex,ey]=P(R.dens[0],at(ch,D.chrlen[ch])); d+=`L${ex} ${ey}Z`;
    el(s,'path',{d,fill:c.data2,opacity:.34});
  });
  CHR22.forEach(ch=>{const [a0,a1]=arcs[ch]; arc(R.la[0],R.la[1],a0,a1,a.n,.5)});
  D.la_segments.forEach(g=>{
    const col=(_P1&&g.anc===_P1)?a.s:(g.anc===_P0?(a.n||c.bg):c.fd);
    const a0=at(g.chrom,g.start),a1=at(g.chrom,g.end);
    const h=(R.la[1]-R.la[0])/2, r0=g.hap===1?R.la[0]+h:R.la[0];
    const p=arc(r0,r0+h,a0,Math.max(a1,a0+0.0016),col,.95);
    tip(p,`chr${g.chrom}:${fmt(g.start)}-${fmt(g.end)} · ${g.anc} · hap ${g.hap}`);
  });
  CHR22.forEach(ch=>{const [a0,a1]=arcs[ch]; arc(R.arch[0],R.arch[1],a0,a1,c.track)});
  D.archaic.forEach(g=>{
    const col=g.src==='Denisovan'?c.data2:a.a;
    const p=arc(R.arch[0],R.arch[1],at(g.chrom,g.start),Math.max(at(g.chrom,g.end),at(g.chrom,g.start)+0.0016),col,g.hom?1:.7);
    tip(p,`chr${g.chrom}:${fmt(g.start)}-${fmt(g.end)} · ${g.src} · ${g.mb} Mb${g.hom?' · homozygous':''}`);
  });
  CHR22.forEach(ch=>{const [a0,a1]=arcs[ch]; arc(R.roh[0],R.roh[1],a0,a1,c.track)});
  D.roh.forEach(r=>{ if(!arcs[r.chrom]) return;
    const p=arc(R.roh[0],R.roh[1],at(r.chrom,r.start),Math.max(at(r.chrom,r.end),at(r.chrom,r.start)+0.0016),c.ink,r.mb>=1.5?.95:.45);
    tip(p,`chr${r.chrom}:${fmt(r.start)}-${fmt(r.end)} · ${r.mb} Mb`);
  });
  const DEL=(D.sv_gene_dels||[]).filter(g=>CHR22.indexOf(String(g.chrom))>=0).map(g=>[g.chrom,g.pos,g.gene,g.frac]);
  DEL.forEach(([ch,pos,name,drop])=>{
    const t=at(ch,pos); const [x0,y0]=P(R.roh[0]-3,t),[x1,y1]=P(R.roh[0]-drop,t);
    el(s,'line',{x1:x0,y1:y0,x2:x1,y2:y1,stroke:c.faint,'stroke-width':.8});
    const right=Math.cos(t)>=0;
    txt(s,{x:x1+(right?4:-4),y:y1,'text-anchor':right?'start':'end','dominant-baseline':'central','font-size':9,'font-weight':700,'font-style':'italic',fill:c.lab,stroke:c.bg,'stroke-width':2.6,'paint-order':'stroke'},name);
  });
  txt(s,{x:CX,y:CY-16,'text-anchor':'middle','font-size':12,'font-weight':700,fill:c.muted,'letter-spacing':'.16em'},NAME().toUpperCase());
  txt(s,{x:CX,y:CY+14,'text-anchor':'middle','font-size':30,'font-weight':800,fill:c.ink,...NUM},(((D.chrlen?Object.entries(D.chrlen).filter(([c])=>/^\d+$/.test(c)&&+c<23).reduce((a,e)=>a+(e[1]||0),0):0)/1e9).toFixed(2))+' Gb');
  txt(s,{x:CX,y:CY+34,'text-anchor':'middle','font-size':10,fill:c.muted},zh()?'22 条常染色体':'22 autosomes');
  const RL=[[R.dens,zh()?'变异密度':'variant density',c.data2,.34],
            [R.la,zh()?'局部祖源':'local ancestry',a.n,.5],
            [R.arch,zh()?'古老人类片段':'archaic segments',a.a,1],
            [R.roh,zh()?'纯合区':'homozygous runs',c.ink,.95]];
  const kx=6, ky=18;
  RL.forEach(([rr,lab,col,op],i)=>{const y=ky+i*18;
    el(s,'rect',{x:kx,y:y-6,width:24,height:9,rx:1.5,fill:c.track});
    el(s,'rect',{x:kx+(i===0?0:5),y:y-6,width:i===0?24:11,height:9,rx:1.5,fill:col,opacity:op});
    txt(s,{x:kx+30,y:y+2,'font-size':9.5,fill:c.lab},lab)});
  const yo=ky+RL.length*18;
  el(s,'rect',{x:kx,y:yo-6,width:24,height:9,rx:1.5,fill:a.n,opacity:.5});
  el(s,'rect',{x:kx+5,y:yo-6,width:11,height:9,rx:1.5,fill:a.s});
  txt(s,{x:kx+30,y:yo+2,'font-size':9.5,fill:c.lab},zh()?'南方东亚片段':'southern segment');
});

/* ── archaic gene families ── */
reveal('afam',(s,c)=>{
  const a=AC(), rows=D.archaic_families.slice(0,10), x0=240,x1=760,top=16,rh=25;
  const X=v=>x0+(x1-x0)*v/100;
  [0,25,50,75,100].forEach(v=>{el(s,'line',{x1:X(v),y1:top-8,x2:X(v),y2:top+rows.length*rh-6,stroke:c.grid,'stroke-width':.6});
    txt(s,{x:X(v),y:top+rows.length*rh+8,'text-anchor':'middle','font-size':8,fill:c.faint},v+'%')});
  txt(s,{x:x0,y:8,'font-size':8,'font-weight':600,fill:c.faint,'letter-spacing':'.08em'},zh()?'该基因家族有多少比例的基因落在渗入片段内':'SHARE OF THE FAMILY\'S GENES THAT FALL INSIDE AN INTROGRESSED SEGMENT');
  rows.forEach((r,i)=>{const y=top+i*rh+6, w=X(r.pct)-x0, hi=r.pct>=25;
    txt(s,{x:x0-10,y:y+3,'text-anchor':'end','font-size':10,'font-weight':hi?800:600,fill:hi?c.ink:c.lab},zh()?r.family:(r.family_en||r.family));
    const b=el(s,'rect',{x:x0,y:y-6,width:Math.max(w,1.5),height:12,rx:2,fill:hi?a.a:c.data2,opacity:hi?1:.7,class:'fade',style:`animation-delay:${i*.05}s`});
    tip(b,`${zh()?r.family:(r.family_en||r.family)}: ${r.carried}/${r.total} · ${r.examples}`);
    txt(s,{x:X(r.pct)+7,y:y+3,'font-size':9.5,'font-weight':800,fill:c.ink},`${r.carried}/${r.total}`)});
  foot(s,c,900,296,zh()?`片段总共覆盖 ${D.archaic_genes.n_genes} 个蛋白编码基因，随机打乱片段位置的期望是 ${Math.round(D.archaic_genes.null_mean)} 个：总数不高于随机，特定家族才是重点`
                     :`the segments cover ${D.archaic_genes.n_genes} protein-coding genes against ${Math.round(D.archaic_genes.null_mean)} expected by chance: the total is unremarkable, the families are the point`);
});

/* ── behaviour PRS ── */
reveal('behaviour',(s,c)=>{
  const a=AC(), rows=[...D.behaviour].sort((p,q)=>q.pct_EAS-p.pct_EAS), x0=230,x1=800,top=24,rh=23;
  const X=v=>x0+(x1-x0)*v/100;
  const H=top+rows.length*rh+44; s.setAttribute('viewBox','0 0 900 '+H);
  el(s,'rect',{x:X(25),y:top-10,width:X(75)-X(25),height:rows.length*rh+2,fill:c.track,opacity:.6});
  [0,25,50,75,100].forEach(v=>{el(s,'line',{x1:X(v),y1:top-10,x2:X(v),y2:top+rows.length*rh-4,stroke:v===50?c.faint:c.grid,'stroke-width':v===50?.9:.5});
    txt(s,{x:X(v),y:top+rows.length*rh+10,'text-anchor':'middle','font-size':8,fill:c.faint},v)});
  txt(s,{x:x0,y:10,'font-size':8,'font-weight':600,fill:c.faint,'letter-spacing':'.08em'},zh()?`在 ${(D.pop||{}).n_super||'—'} 名${_spZh()}参考个体中的百分位`:`PERCENTILE AMONG THE ${(D.pop||{}).n_super||'—'} ${(((D.pop||{}).superpop)||'EAS')} REFERENCE INDIVIDUALS`);
  rows.forEach((r,i)=>{const y=top+i*rh+6, eas=r.panel==='EAS';
    txt(s,{x:x0-12,y:y+3,'text-anchor':'end','font-size':10,'font-weight':eas?700:500,fill:eas?c.ink:c.muted},zh()?r.zh:r.en);
    const mx=X(r.pct_EAS);
    if(eas){const d=el(s,'polygon',{points:`${mx},${y-6.4} ${mx+6.4},${y} ${mx},${y+6.4} ${mx-6.4},${y}`,fill:a.n,stroke:c.bg,'stroke-width':1.2,class:'pop'});
      d.style.animationDelay=(i*25)+'ms'; tip(d,`${r.en}: ${r.pct_EAS} pct · z ${r.z} · coverage ${r.coverage}%`)}
    else {const d=el(s,'circle',{cx:mx,cy:y,r:4.5,fill:c.bg,stroke:a.s,'stroke-width':1.8,class:'pop'});
      d.style.animationDelay=(i*25)+'ms'; tip(d,`${r.en}: ${r.pct_EAS} pct · z ${r.z} · coverage ${r.coverage}% · European-derived`)}
    txt(s,{x:x1+12,y:y+3,'font-size':9.5,'font-weight':800,fill:eas?c.ink:c.muted},r.pct_EAS.toFixed(0));
    txt(s,{x:x1+36,y:y+3,'font-size':8,fill:c.faint},r.coverage+'%')});
  txt(s,{x:x1+12,y:top-10,'font-size':7.5,'font-weight':600,fill:c.faint,'letter-spacing':'.06em'},zh()?'百分位 · 覆盖':'PCT · COV');
  // 两种来源的分数**都用东亚参考分布**定位（26 的输出字段就叫 pct_EAS）。欧洲训练的评分在这里只是
  // 迁移性对照：它的百分位不代表该个体在欧洲人群内的位置，图注必须说出来，否则读者会那样理解。
  foot(s,c,900,H-8,zh()?`菱形 = 东亚人群训练 · 空心圆 = 欧洲人群训练（其百分位同样在东亚参考分布下计算，不代表在欧洲人群内的位置；跨人群迁移性差，只作对照）· 灰带 = 中间一半的东亚参考个体`
                        :`diamond = trained in East Asians · hollow circle = trained in Europeans (its percentile is also computed against the East Asian reference, not a position within Europeans; cross-population transfer is poor, shown for contrast only) · band = middle half of the East Asian reference`);
});

renderAll();
