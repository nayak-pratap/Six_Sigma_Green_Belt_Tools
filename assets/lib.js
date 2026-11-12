/*MATH-START*/
function erfcc(x) {
  var z = Math.abs(x);
  var t = 1 / (1 + 0.5 * z);
  var r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
    t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
  return x >= 0 ? r : 2 - r;
}
function normCdf(x) { return 0.5 * erfcc(-x / Math.SQRT2); }
function normInv(p) {
  if (p <= 0) return -Infinity;
  if (p >= 1) return Infinity;
  var a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
  var b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01];
  var c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
  var d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];
  var pl = 0.02425, q, r, x;
  if (p < pl) {
    q = Math.sqrt(-2 * Math.log(p));
    x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  } else if (p <= 1 - pl) {
    q = p - 0.5; r = q * q;
    x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
  } else {
    q = Math.sqrt(-2 * Math.log(1 - p));
    x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  }
  // one Halley refinement step
  var e = normCdf(x) - p;
  var u = e * Math.sqrt(2 * Math.PI) * Math.exp(x * x / 2);
  return x - u / (1 + x * u / 2);
}
function gammaln(x) {
  var cof = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
  var y = x, tmp = x + 5.5;
  tmp -= (x + 0.5) * Math.log(tmp);
  var ser = 1.000000000190015;
  for (var j = 0; j < 6; j++) { y += 1; ser += cof[j] / y; }
  return -tmp + Math.log(2.5066282746310005 * ser / x);
}
function betacf(a, b, x) {
  var MAXIT = 300, EPS = 3e-12, FPMIN = 1e-300;
  var qab = a + b, qap = a + 1, qam = a - 1, c = 1, d = 1 - qab * x / qap;
  if (Math.abs(d) < FPMIN) d = FPMIN;
  d = 1 / d;
  var h = d;
  for (var m = 1; m <= MAXIT; m++) {
    var m2 = 2 * m;
    var aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d; h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    var del = d * c; h *= del;
    if (Math.abs(del - 1) < EPS) break;
  }
  return h;
}
function betai(a, b, x) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  var bt = Math.exp(gammaln(a + b) - gammaln(a) - gammaln(b) + a * Math.log(x) + b * Math.log(1 - x));
  if (x < (a + 1) / (a + b + 2)) return bt * betacf(a, b, x) / a;
  return 1 - bt * betacf(b, a, 1 - x) / b;
}
function tTwoSidedP(t, df) { return betai(df / 2, 0.5, df / (df + t * t)); }
function tCdf(t, df) {
  var p2 = tTwoSidedP(t, df);
  return t >= 0 ? 1 - p2 / 2 : p2 / 2;
}
function tInv(p, df) {
  var lo = -200, hi = 200;
  for (var i = 0; i < 200; i++) {
    var mid = (lo + hi) / 2;
    if (tCdf(mid, df) < p) lo = mid; else hi = mid;
  }
  return (lo + hi) / 2;
}

function sigmaLevel(defects, units, opps, shift) {
  var dpu = defects / units;
  var dpmo = defects / (units * opps) * 1e6;
  var yld = 1 - defects / (units * opps);
  var z;
  if (dpmo <= 0) z = Infinity;
  else if (dpmo >= 1e6) z = -Infinity;
  else z = normInv(1 - dpmo / 1e6);
  return { dpu: dpu, dpmo: dpmo, yield: yld, sigma: z + (shift ? 1.5 : 0), z: z };
}
function capability(usl, lsl, mean, sd) {
  var cp = (usl - lsl) / (6 * sd);
  var cpu = (usl - mean) / (3 * sd);
  var cpl = (mean - lsl) / (3 * sd);
  var cpk = Math.min(cpu, cpl);
  var ppmLow = normCdf((lsl - mean) / sd) * 1e6;
  var ppmHigh = normCdf(-(usl - mean) / sd) * 1e6;
  return { cp: cp, cpu: cpu, cpl: cpl, cpk: cpk, ppmLow: ppmLow, ppmHigh: ppmHigh, ppm: ppmLow + ppmHigh };
}
var CC = { 2: [1.880, 0, 3.267, 1.128], 3: [1.023, 0, 2.574, 1.693], 4: [0.729, 0, 2.282, 2.059], 5: [0.577, 0, 2.114, 2.326],
  6: [0.483, 0, 2.004, 2.534], 7: [0.419, 0.076, 1.924, 2.704], 8: [0.373, 0.136, 1.864, 2.847], 9: [0.337, 0.184, 1.816, 2.970], 10: [0.308, 0.223, 1.777, 3.078] };
function xbarR(n, grand, rbar) {
  var k = CC[n];
  return { xUcl: grand + k[0] * rbar, xCl: grand, xLcl: grand - k[0] * rbar, rUcl: k[2] * rbar, rCl: rbar, rLcl: k[1] * rbar, sigma: rbar / k[3] };
}
function imr(mean, mrbar) {
  return { xUcl: mean + 2.66 * mrbar, xCl: mean, xLcl: mean - 2.66 * mrbar, rUcl: 3.267 * mrbar, rCl: mrbar, rLcl: 0, sigma: mrbar / 1.128 };
}
function pChart(pbar, n) {
  var s = Math.sqrt(pbar * (1 - pbar) / n);
  return { ucl: pbar + 3 * s, cl: pbar, lcl: Math.max(0, pbar - 3 * s) };
}
function npChart(pbar, n) {
  var s = Math.sqrt(n * pbar * (1 - pbar));
  return { ucl: n * pbar + 3 * s, cl: n * pbar, lcl: Math.max(0, n * pbar - 3 * s) };
}
function cChart(cbar) {
  var s = Math.sqrt(cbar);
  return { ucl: cbar + 3 * s, cl: cbar, lcl: Math.max(0, cbar - 3 * s) };
}
function welch(m1, s1, n1, m2, s2, n2) {
  var a = s1 * s1 / n1, b = s2 * s2 / n2;
  var se = Math.sqrt(a + b);
  var t = (m1 - m2) / se;
  var df = (a + b) * (a + b) / (a * a / (n1 - 1) + b * b / (n2 - 1));
  var p = tTwoSidedP(t, df);
  var tc = tInv(0.975, df);
  var diff = m1 - m2;
  return { t: t, df: df, p: p, diff: diff, se: se, lo: diff - tc * se, hi: diff + tc * se };
}
function grr(ev, av, pv) {
  var g = Math.sqrt(ev * ev + av * av);
  var tv = Math.sqrt(g * g + pv * pv);
  return { grr: g, tv: tv, pctGrr: g / tv * 100, pctEv: ev / tv * 100, pctAv: av / tv * 100, pctPv: pv / tv * 100, ndc: Math.floor(1.41 * pv / g) };
}
function planSchedule(weeks, hoursPerWeek) {
  var weights = [['Foundations', 15], ['Define', 10], ['Measure', 20], ['Analyze', 25], ['Improve', 15], ['Control', 15]];
  var reviewWeeks = weeks >= 8 ? 2 : 1;
  var learnWeeks = weeks - reviewWeeks;
  var pool = learnWeeks * hoursPerWeek;
  var queue = weights.map(function (w) { return { name: w[0], left: pool * w[1] / 100 }; });
  var out = [], qi = 0;
  for (var w = 1; w <= learnWeeks; w++) {
    var cap = hoursPerWeek, items = [];
    while (cap > 1e-9 && qi < queue.length) {
      var take = Math.min(cap, queue[qi].left);
      if (take > 1e-9) items.push({ name: queue[qi].name, h: take });
      queue[qi].left -= take; cap -= take;
      if (queue[qi].left <= 1e-9) qi++;
    }
    out.push({ week: w, items: items, kind: 'learn' });
  }
  for (var r = 0; r < reviewWeeks; r++) {
    out.push({ week: learnWeeks + r + 1, items: [{ name: r === reviewWeeks - 1 ? 'Full mock exam and weak-area review' : 'Mixed practice by phase and calculation drills', h: hoursPerWeek }], kind: 'review' });
  }
  return out;
}
/*MATH-END*/

/* ---------- helpers shared by every page ---------- */
function $(s, r) { return (r || document).querySelector(s); }
function $$(s, r) { return Array.prototype.slice.call((r || document).querySelectorAll(s)); }
function esc(s) { return String(s).replace(/[&<>"]/g, function (c) { return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]; }); }
function fmt(x, d) { if (x === Infinity || x === -Infinity || x !== x) return '—'; return Number(x).toFixed(d === undefined ? 2 : d); }
function val(id) { var e = document.getElementById(id); return e ? parseFloat(e.value) : NaN; }
function kpi(label, v, cls) { return '<div class="kpi ' + (cls || '') + '"><span>' + esc(label) + '</span><b>' + v + '</b></div>'; }
function parseNums(s) { return String(s).split(/[\s,;]+/).map(parseFloat).filter(function (x) { return x === x && isFinite(x); }); }
function mean(a) { return a.reduce(function (x, y) { return x + y; }, 0) / a.length; }
function sdev(a) { var m = mean(a); return Math.sqrt(a.reduce(function (x, y) { return x + (y - m) * (y - m); }, 0) / (a.length - 1)); }
function quant(sorted, p) { var h = (sorted.length - 1) * p, lo = Math.floor(h), hi = Math.ceil(h); return sorted[lo] + (sorted[hi] - sorted[lo]) * (h - lo); }
function niceTicks(lo, hi, n) {
  n = n || 5; var span = hi - lo; if (!(span > 0)) return [lo];
  var raw = span / n, mag = Math.pow(10, Math.floor(Math.log(raw) / Math.LN10)), f = raw / mag, step = (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * mag, out = [];
  for (var v = Math.ceil(lo / step - 1e-9) * step; v <= hi + step * 1e-9; v += step) out.push(+v.toPrecision(12));
  return out;
}
function frame(W, H, x0, x1, y0, y1, o) {
  o = o || {};
  var ml = o.ml || 52, mr = o.mr || 14, mt = o.mt || 14, mb = o.mb || 40;
  var f = { W: W, H: H, ml: ml, mr: mr, mt: mt, mb: mb,
    xs: function (v) { return ml + (v - x0) / ((x1 - x0) || 1) * (W - ml - mr); },
    ys: function (v) { return H - mb - (v - y0) / ((y1 - y0) || 1) * (H - mb - mt); } };
  var t = 'style="font:11px var(--font-mono);fill:var(--ink-2)"', s = '';
  niceTicks(y0, y1, 5).forEach(function (v) {
    s += '<line x1="' + ml + '" x2="' + (W - mr) + '" y1="' + f.ys(v) + '" y2="' + f.ys(v) + '" style="stroke:var(--line)"/><text x="' + (ml - 8) + '" y="' + (f.ys(v) + 4) + '" text-anchor="end" ' + t + '>' + (+v.toFixed(6)) + '</text>';
  });
  if (o.xticks !== false) niceTicks(x0, x1, o.xn || 6).forEach(function (v) { s += '<text x="' + f.xs(v) + '" y="' + (H - mb + 16) + '" text-anchor="middle" ' + t + '>' + (+v.toFixed(6)) + '</text>'; });
  if (o.xl) s += '<text x="' + ((ml + W - mr) / 2) + '" y="' + (H - 4) + '" text-anchor="middle" ' + t + '>' + esc(o.xl) + '</text>';
  if (o.yl) s += '<text transform="translate(12 ' + ((mt + H - mb) / 2) + ') rotate(-90)" text-anchor="middle" ' + t + '>' + esc(o.yl) + '</text>';
  f.base = s; f.t = t;
  return f;
}
function svgWrap(W, H, inner, label) { return '<svg viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="' + esc(label) + '" style="width:100%;height:auto">' + inner + '</svg>'; }

/* Wire a page: fill example values, recompute on any input, and make Reset restore them. */
function setupTool(defaults, update, customReset) {
  function apply() {
    if (customReset) { customReset(); return; }
    Object.keys(defaults).forEach(function (id) { var e = document.getElementById(id); if (e) e.value = defaults[id]; });
  }
  apply();
  var host = document.getElementById('tool');
  host.addEventListener('input', function () { update(); });
  host.addEventListener('change', function () { update(); });
  document.getElementById('resetBtn').addEventListener('click', function () { apply(); update(); });
  update();
}
function copyText(text, btn) {
  var old = btn.textContent;
  function done(m) { btn.textContent = m; setTimeout(function () { btn.textContent = old; }, 1800); }
  try { navigator.clipboard.writeText(text).then(function () { done('Copied'); }, function () { done('Copy blocked: select the text'); }); }
  catch (e) { done('Copy blocked: select the text'); }
}

/* Decision-tree chooser used by the test and chart pages */
function makeChooser(cfg, kind) {
  var st = { a: cfg.first.options[0][0], b: 0 };
  var el = document.createElement('div'); el.className = 'chooser';
  function draw() {
    var h = '<h3>' + esc(cfg.title) + '</h3><div><div class="q-label">' + esc(cfg.first.label) + '</div><div class="opts">';
    cfg.first.options.forEach(function (o) { h += '<button class="opt" type="button" data-l="a" data-v="' + o[0] + '" aria-pressed="' + (st.a === o[0]) + '">' + esc(o[1]) + '</button>'; });
    h += '</div></div><div><div class="q-label">What are you trying to do?</div><div class="opts">';
    cfg.results[st.a].forEach(function (r, i) { h += '<button class="opt" type="button" data-l="b" data-v="' + i + '" aria-pressed="' + (st.b === i) + '">' + esc(r[0]) + '</button>'; });
    h += '</div></div>';
    var r = cfg.results[st.a][st.b];
    if (kind === 'test') h += '<div class="result"><b class="t">' + esc(r[1]) + '</b><dl><dt>If not normal</dt><dd>' + esc(r[2]) + '</dd><dt>Check</dt><dd>' + esc(r[3]) + '</dd><dt>Example</dt><dd>' + esc(r[4]) + '</dd></dl></div>';
    else h += '<div class="result"><b class="t">' + esc(r[1]) + '</b><dl><dt>Use when</dt><dd>' + esc(r[2]) + '</dd><dt>Limits</dt><dd>' + esc(r[3]) + '</dd></dl></div>';
    el.innerHTML = h;
  }
  el.addEventListener('click', function (e) {
    var b = e.target.closest('.opt'); if (!b) return;
    if (b.getAttribute('data-l') === 'a') { st.a = b.getAttribute('data-v'); st.b = 0; } else { st.b = parseInt(b.getAttribute('data-v'), 10); }
    draw();
  });
  el.reset = function () { st.a = cfg.first.options[0][0]; st.b = 0; draw(); };
  draw();
  return el;
}

/* theme toggle, shares the key used by the framework page */
(function () {
  try { var th = localStorage.getItem('gb-theme'); if (th) document.documentElement.setAttribute('data-theme', th); } catch (e) { }
  document.addEventListener('DOMContentLoaded', function () {
    var b = document.getElementById('themeBtn'); if (!b) return;
    b.addEventListener('click', function () {
      var cur = document.documentElement.getAttribute('data-theme') || ((window.matchMedia && matchMedia('(prefers-color-scheme: dark)').matches) ? 'dark' : 'light');
      var n = cur === 'dark' ? 'light' : 'dark';
      document.documentElement.setAttribute('data-theme', n);
      try { localStorage.setItem('gb-theme', n); } catch (e) { }
    });
  });
})();

/* ---------- generic engines used by the template-style tools ---------- */
function nz(v) { v = parseFloat(v); return isFinite(v) ? v : 0; }
function hbars(items, fmtv) {
  var mx = Math.max.apply(null, items.map(function (i) { return i.v; }).concat([1e-9]));
  return '<div class="meters" style="display:flex;flex-direction:column;gap:6px">' + items.map(function (i) {
    return '<div class="prow"><span>' + esc(i.l) + '</span><div class="bar"><i style="width:' + Math.round(i.v / mx * 100) + '%' + (i.c ? ';background:' + i.c : '') + '"></i></div><span class="n">' + (fmtv ? fmtv(i.v) : fmt(i.v, 1)) + '</span></div>';
  }).join('') + '</div>';
}
function gridTable(cfg) {
  var host = document.getElementById(cfg.host), rows = [], extra = cfg.extra || [];
  function cell(c, v) {
    if (c.t === 'sel') return '<select data-k="' + c.k + '" aria-label="' + esc(c.l) + '">' + c.o.map(function (o) {
      var ov = Array.isArray(o) ? o[0] : o, ol = Array.isArray(o) ? o[1] : o;
      return '<option value="' + esc(ov) + '"' + (String(v) === String(ov) ? ' selected' : '') + '>' + esc(ol) + '</option>';
    }).join('') + '</select>';
    if (c.t === 'num') return '<input class="fmea-in" type="number" step="any"' + (c.min !== undefined ? ' min="' + c.min + '"' : '') + (c.max !== undefined ? ' max="' + c.max + '"' : '') + ' data-k="' + c.k + '" value="' + (v === undefined || v === null ? '' : v) + '" aria-label="' + esc(c.l) + '">';
    return '<input class="fmea-name" type="text" data-k="' + c.k + '" value="' + esc(v === undefined ? '' : v) + '" aria-label="' + esc(c.l) + '" style="min-width:' + (c.w || 120) + 'px">';
  }
  function draw() {
    var h = '<div class="tbl-wrap"><table><thead><tr>' + cfg.cols.map(function (c) { return '<th>' + esc(c.l) + '</th>'; }).join('') + extra.map(function (e) { return '<th>' + esc(e) + '</th>'; }).join('') + '<th></th></tr></thead><tbody>';
    rows.forEach(function (r, i) {
      h += '<tr data-i="' + i + '">' + cfg.cols.map(function (c) { return '<td>' + cell(c, r[c.k]) + '</td>'; }).join('') + extra.map(function (e, j) { return '<td data-x="' + j + '"></td>'; }).join('') + '<td><button class="del" type="button" data-del="' + i + '" aria-label="Delete row">×</button></td></tr>';
    });
    if (!rows.length) h += '<tr><td colspan="' + (cfg.cols.length + extra.length + 1) + '" class="muted">No rows. Add one.</td></tr>';
    h += '</tbody></table></div><div class="acts"><button class="btn small" type="button" data-add>' + esc(cfg.addLabel || 'Add row') + '</button><button class="btn small" type="button" data-clear>Clear all rows</button><button class="btn small" type="button" data-copy>Copy as text</button></div>';
    host.innerHTML = h; recompute();
  }
  function recompute() {
    var res = cfg.compute ? cfg.compute(rows) : {};
    (res.cells || []).forEach(function (cs, i) {
      var tr = host.querySelector('tr[data-i="' + i + '"]'); if (!tr) return;
      cs.forEach(function (html, j) { var td = tr.querySelector('td[data-x="' + j + '"]'); if (td) td.innerHTML = html; });
    });
    var out = document.getElementById(cfg.out || 'out'); if (out) out.innerHTML = res.out || '';
  }
  function onEdit(e) {
    var el = e.target.closest('[data-k]'); if (!el) return;
    var i = +el.closest('tr').getAttribute('data-i'), k = el.getAttribute('data-k'), c = cfg.cols.filter(function (x) { return x.k === k; })[0], v = el.value;
    if (c.t === 'num') { v = v === '' ? '' : parseFloat(v); if (v !== '' && !isFinite(v)) v = ''; }
    rows[i][k] = v; recompute();
  }
  host.addEventListener('input', onEdit); host.addEventListener('change', onEdit);
  host.addEventListener('click', function (e) {
    var t = e.target;
    if (t.hasAttribute('data-del')) { rows.splice(+t.getAttribute('data-del'), 1); draw(); }
    else if (t.hasAttribute('data-add')) { var nr = {}; cfg.cols.forEach(function (c) { nr[c.k] = c.d !== undefined ? c.d : (c.t === 'num' ? '' : (c.t === 'sel' ? (Array.isArray(c.o[0]) ? c.o[0][0] : c.o[0]) : '')); }); rows.push(nr); draw(); }
    else if (t.hasAttribute('data-clear')) { rows = []; draw(); }
    else if (t.hasAttribute('data-copy')) {
      copyText(rows.map(function (r) { return cfg.cols.map(function (c) { return c.l + ': ' + r[c.k]; }).join(' | '); }).join('\n'), t);
    }
  });
  return {
    reset: function () { rows = cfg.example.map(function (r) { var o = {}; Object.keys(r).forEach(function (k) { o[k] = r[k]; }); return o; }); draw(); },
    recompute: recompute, rows: function () { return rows; }
  };
}
/* scored audit: groups of statements, one select each */
function auditTool(cfg) {
  var host = document.getElementById(cfg.host), S = [];
  var mx = Math.max.apply(null, cfg.scale.map(function (s) { return s[0]; }));
  function draw() {
    host.innerHTML = cfg.groups.map(function (g, gi) {
      return '<div style="margin-bottom:18px"><h3 style="margin-bottom:8px">' + esc(g.n) + '</h3><div class="tbl-wrap"><table><tbody>' + g.items.map(function (it, ii) {
        return '<tr><td>' + esc(it) + '</td><td style="width:210px"><select data-g="' + gi + '" data-i="' + ii + '" aria-label="' + esc(it) + '">' + cfg.scale.map(function (s) { return '<option value="' + s[0] + '"' + (S[gi][ii] === s[0] ? ' selected' : '') + '>' + esc(s[1]) + '</option>'; }).join('') + '</select></td></tr>';
      }).join('') + '</tbody></table></div></div>';
    }).join('');
    compute();
  }
  function compute() {
    var tot = 0, max = 0, per = cfg.groups.map(function (g, gi) { var s = S[gi].reduce(function (a, b) { return a + b; }, 0), m = g.items.length * mx; tot += s; max += m; return { l: g.n, v: m ? s / m * 100 : 0 }; });
    var pct = max ? tot / max * 100 : 0, band = cfg.bands.filter(function (b) { return pct >= b[0]; })[0] || cfg.bands[cfg.bands.length - 1];
    var low = per.slice().sort(function (a, b) { return a.v - b.v; })[0];
    document.getElementById('out').innerHTML = '<div class="kpis">' + kpi('Overall score', fmt(pct, 0) + '%', band[2]) + kpi('Lowest area', esc(low.l), 'warn') + '</div>' +
      hbars(per.map(function (p) { return { l: p.l, v: p.v }; }), function (v) { return fmt(v, 0) + '%'; }) +
      '<div class="verdict ' + band[2] + '" style="margin-top:12px"><b>' + esc(band[1]) + '</b> ' + esc(cfg.advice(low.l)) + '</div>';
  }
  host.addEventListener('change', function (e) { var s = e.target.closest('select[data-g]'); if (!s) return; S[+s.getAttribute('data-g')][+s.getAttribute('data-i')] = parseFloat(s.value); compute(); });
  return { reset: function () { S = cfg.example.map(function (r) { return r.slice(); }); draw(); } };
}
