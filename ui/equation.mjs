// Equations split into their physical terms.
//
// An equation is shown as a row of term cards. Each card holds one term's maths,
// its name, its live numerical value with units, and a colour that matches the
// curve or element it produces in the plots. Hovering a term highlights its curve.
// A "where" block defines every symbol with its value; a derivation lists steps.
//
//   const eq = createEquation({
//     lhs: { mathml: '<mi>V</mi><mo>(</mo><mi>r</mi><mo>)</mo>', label: 'potential', value: () => V, unit: 'GeV' },
//     terms: [
//       { id: 'coul', mathml: '<mo>−</mo><mfrac>…</mfrac>', label: 'Coulomb (gluon exchange)', color: '--body-1', value: () => Vc, unit: 'GeV' },
//       { id: 'lin', sign: '+', mathml: '<mi>σ</mi><mi>r</mi>', label: 'string (confinement)', color: '--body-0', value: () => s * r, unit: 'GeV' },
//     ],
//     where: [{ mathml: '<mi>σ</mi>', meaning: 'string tension', value: () => 0.926, unit: 'GeV/fm' }],
//     link: plot,
//   });
//   eq.update();   // refresh live values (call every UI tick)

import { el, sci } from './instrument.mjs';

const m = inner => `<math display="inline">${inner}</math>`;

function fmt(v, digits = 4) {
  if (typeof v === 'string') return v;
  if (v == null || Number.isNaN(v)) return '—';
  const a = Math.abs(v);
  if (a !== 0 && (a < 1e-3 || a >= 1e5)) return sci(v, digits - 1);
  return (+v.toPrecision(digits)).toString().replace('-', '−');
}

export function createEquation({ title, lhs, op = '=', terms = [], where = [], note, link, onFocus } = {}) {
  const valueEls = [];
  const focusTerm = id => {
    for (const c of root.querySelectorAll('.psim-term')) c.classList.toggle('dim', id != null && c.dataset.id !== id);
    const targets = Array.isArray(link) ? link : link ? [link] : [];
    targets.forEach(p => p.highlight(id));
    onFocus?.(id);
  };

  const card = (t, kind) => {
    const v = el('span', { class: 'v' });
    if (t.value) valueEls.push([v, t]);
    const c = el('div', {
      class: `psim-term ${kind}`, 'data-id': t.id ?? '', tabindex: t.id ? 0 : null,
      onmouseenter: () => t.id && focusTerm(t.id), onmouseleave: () => t.id && focusTerm(null),
      onfocus: () => t.id && focusTerm(t.id), onblur: () => t.id && focusTerm(null),
    },
    el('div', { class: 'mx', html: m(t.mathml) }),
    el('i', { class: 'bar', style: t.color ? `background:var(${t.color})` : '' }),
    t.label ? el('span', { class: 'lb' }, t.label) : null,
    v);
    return c;
  };

  const row = el('div', { class: 'psim-eqrow' });
  if (lhs) row.append(card(lhs, 'lhs'), el('span', { class: 'op' }, op));
  terms.forEach((t, i) => {
    if (i > 0) row.append(el('span', { class: 'op' }, t.sign ?? '+'));
    row.append(card(t, 'rhs'));
  });

  const whereEl = where.length ? el('div', { class: 'psim-where' },
    el('span', { class: 'kw' }, 'where'),
    ...where.map(w => {
      const v = el('span', { class: 'v' });
      if (w.value) valueEls.push([v, w]);
      return el('div', { class: 'w' }, el('span', { class: 'mx', html: m(w.mathml) }), el('span', { class: 'mean' }, w.meaning), v);
    })) : null;

  const root = el('section', { class: 'psim-equation' },
    title ? el('h3', {}, title) : null, row, whereEl, note ? el('p', { class: 'psim-hint' }, note) : null);

  const api = {
    el: root,
    update() {
      for (const [node, t] of valueEls) {
        let v;
        try { v = t.value(); } catch { v = NaN; }
        node.textContent = `${fmt(v, t.digits)}${t.unit ? ` ${t.unit}` : ''}`;
      }
    },
    focus: focusTerm,
  };
  api.update();
  return api;
}

/** Numbered derivation: each step is one line of maths plus the reason for it. */
export function createDerivation(steps, { title } = {}) {
  return el('section', { class: 'psim-derivation' },
    title ? el('h3', {}, title) : null,
    el('ol', {}, steps.map(s => el('li', {},
      el('div', { class: 'mx', html: `<math display="inline">${s.mathml}</math>` }),
      s.note ? el('span', { class: 'why' }, s.note) : null))));
}

/** Helpers for writing MathML compactly. */
export const M = {
  i: s => `<mi>${s}</mi>`,
  n: s => `<mn>${s}</mn>`,
  o: s => `<mo>${s}</mo>`,
  // msub/msup take exactly two children: wrap both, so a multi-part exponent such as −x survives
  sub: (b, s) => `<msub><mrow>${b}</mrow><mrow>${s}</mrow></msub>`,
  sup: (b, s) => `<msup><mrow>${b}</mrow><mrow>${s}</mrow></msup>`,
  frac: (a, b) => `<mfrac><mrow>${a}</mrow><mrow>${b}</mrow></mfrac>`,
  sqrt: a => `<msqrt>${a}</msqrt>`,
  row: (...a) => `<mrow>${a.join('')}</mrow>`,
  paren: a => `<mo>(</mo>${a}<mo>)</mo>`,
  /** Vector symbol: Unicode mathematical bold (MathML Core ignores mathvariant). */
  vec: s => `<mi>${{ r: '𝐫', v: '𝐯', E: '𝐄', B: '𝐁', u: '𝐮', F: '𝐅', p: '𝐩', k: '𝐤', x: '𝐱', a: '𝐚' }[s] ?? s}</mi>`,
};
