"use strict";
/* Reading a provenance edge: the parts that are the same wherever an edge is shown.

   The registry viewer shows the edges that reach one structure; the energy graph walks
   them backwards and shows the same edges again as candidate steps.  Both have to make
   the same decisions — what to call a species, when two edges are one chemistry, how a
   number carries its caveat, what a refusal reads like — and two answers to any of those
   is two different accounts of the same registry.

   Pure functions over the shape `/api/structures/{id}[/routes]` returns, and HTML
   strings rather than nodes because both callers build their tables that way.  Nothing
   here touches the DOM, holds state, or knows which page it is on: a page's own
   rendering (which rows are open, what a click does) stays with the page.

   Exposed as one object; a page destructures the names it wants. */

const mofsbuRoutes = (() => {

/* One escaper for the whole app — chrome.js already has it and it uses no `this`. */
const esc = mofsbuChrome.esc;

const fmtE = (v) => v === null || v === undefined ? "—" : Number(v).toFixed(4);

/* The rungs of the fidelity ladder, by their stored integer (mofsbu._types.Fidelity). */
const FIDELITY = {"-1":"heuristic", "0":"raw", "1":"FF", "2":"ML", "3":"xTB", "4":"DFT"};

/* A refusal from the reference scheme is a paragraph; a list of routes has room for a
   clause. The whole of it goes in the title attribute rather than being thrown away. */
const firstClause = (s) => {
  const t = String(s ?? "").split(/[;.]\s/)[0].trim();
  return t.length > 80 ? t.slice(0, 77) + "…" : t;
};

/* The stored label spells out the identity work: every deprotonation orbit, the spin, the
   identity level it was resolved to.  That is the right thing to store and the wrong thing
   to put in six table cells — it is wider than the panel.  The full string stays one hover
   away, and `renderTerms` in the opened row still shows it in full. */
function shortLabel(full){
  return String(full ?? "")
    .replace(/cfg\([^)]*\)/g, "").replace(/sep\([^)]*\)/g, "")   // the orbit bookkeeping
    .replace(/\s*\?L\d+\s*$/, "")                                 // identity level
    .replace(/\s+s\d+\b/g, "")                                    // spin multiplicity
    .replace(/\s{2,}/g, " ").trim();
}

/* The species on each side of one edge's arrow.  Every id is a link because the caller
   can fetch any structure by id, page or no page; the click is caught by delegation. */
function renderTerms(step){
  const terms = (step && step.terms) || [];
  if(!terms.length)
    return `<span class="muted">no species recorded for this edge</span>`;
  const link = (t) => {
    const name = (t.display_label && t.display_label.trim())
              || (t.label && t.label.trim()) || ("#" + t.structure_id);
    const n = Number(t.stoich) > 1 ? esc(t.stoich) + " " : "";
    return `${n}<a href="#" class="term" data-sid="${esc(t.structure_id)}"`
         + ` title="open structure ${esc(t.structure_id)}">${esc(name)}</a>`;
  };
  const left  = terms.filter(t => t.side === "reagent").map(link).join(" + ");
  const right = terms.filter(t => t.side === "product").map(link).join(" + ");
  return (left  || `<span class="muted">no reagent recorded for this edge</span>`)
       + " → "
       + (right || `<span class="muted">no product recorded for this edge</span>`);
}

/* One row per CHEMISTRY: same source, same thing added, one entry.  The choice vector
   distinguishes two placements of those same pieces, which is a geometric fact and not a
   different reaction — dE does not depend on it — so it is counted inside the group
   rather than splitting it.  Collapsed in the BROWSER: the API's route count is what
   `tests/test_ui_controls.py` pins, and it must not move. */
function groupRoutes(routes){
  const seen = new Map(), out = [];
  // Both sides are in the key: edges that CONSUME one structure share their reagents
  // and can still make different products (two isomers of one adduct).
  const side = (terms, s) => terms.filter(t => t.side === s)
    .map(t => Number(t.structure_id)).sort((a, b) => a - b).join(",");
  for(const r of routes){
    const step = (r.steps && r.steps[0]) || null;
    const terms = (step && step.terms) || [];
    const key = side(terms, "reagent") + ">" + side(terms, "product");
    let g = seen.get(key);
    if(!g){ g = {key, rep:r, members:[], vectors:new Set()}; seen.set(key, g); out.push(g); }
    g.members.push(r);
    if(r.choice_vector_digest) g.vectors.add(r.choice_vector_digest);
  }
  for(const g of out){ g.count = g.members.length; g.placements = g.vectors.size; }
  return out;
}

/* A deprotonation adds nothing — it shifts a proton onto the carrier and the carrier
   leaves.  The `leaving` stoichiometry is how many went. */
function protonsShed(step){
  const n = ((step && step.terms) || [])
    .filter(t => t.role === "leaving")
    .reduce((s, t) => s + Number(t.stoich || 1), 0);
  if(!n) return null;
  return `<span class="caveat">&minus;${esc(n)}&nbsp;H<sup>+</sup></span>`;
}

/* A species cell: short enough to read across, with the stored label one hover away and
   the id still the thing you click. */
function speciesCell(term){
  if(!term) return `<span class="muted">—</span>`;
  const full = term.display_label || ("structure " + term.structure_id);
  const n = Number(term.stoich || 1) > 1 ? esc(term.stoich) + "&times; " : "";
  return `${n}<a href="#" class="term" data-sid="${esc(term.structure_id)}"`
       + ` title="${esc(full)}">${esc(shortLabel(full))}</a>`;
}

/* A dE here is never isodesmic (an assembly step forms a metal-donor bond, D17), so the
   qualification travels with the number rather than being left to the reader to
   remember.  The number and its qualifiers come back apart, because how much width
   they have is the caller's problem: a wide panel sets them on one line, a six-column
   picker stacks them.  What they SAY is decided once, here. */
function energyParts(r){
  const s = (r.steps && r.steps[0]) || {};
  if(!r.can_price)
    return {value: null, why: firstClause(r.why_not) || "no energy", notes: []};
  const notes = [];
  if(s.isodesmic === false)        notes.push("not isodesmic");
  else if(s.isodesmic !== true)    notes.push("isodesmic unchecked");
  if(s.all_converged === false)    notes.push("unconverged");
  return {value: fmtE(r.total_dE), why: null, notes};
}

function routeEnergy(r){
  const p = energyParts(r);
  if(p.value === null) return `<span class="muted">— (${esc(p.why)})</span>`;
  return p.value + p.notes.map(n => ` <span class="caveat">· ${esc(n)}</span>`).join("");
}

function routeRung(r){
  const s = (r.steps && r.steps[0]) || {};
  if(s.fidelity === null || s.fidelity === undefined)
    return `<span class="muted">— (${r.can_price ? "no rung reported" : "not priced"})</span>`;
  return esc(FIDELITY[String(s.fidelity)] || s.fidelity_name || s.fidelity);
}

function routeCaveats(r){
  const s = (r.steps && r.steps[0]) || {};
  const codes = s.caveats || [];
  if(!codes.length)
    return s.isodesmic === true || s.isodesmic === false
      ? `<span class="muted">—</span>`
      : `<span class="muted">— (reference quality not checked)</span>`;
  return `<span class="caveat">${codes.map(c => esc(String(c).replace(/_/g, " ")))
    .join(", ")}</span>`;
}

function routeWhen(r){
  const d = (r.created_at || "").slice(0, 10);
  return d ? esc(d) : `<span class="muted">— (not recorded)</span>`;
}

function routeHead(r){
  return `#${esc(r.id)} · ${esc(r.kind)} · depth ${esc(r.depth ?? "— (M8)")}`
       + (r.note ? " · " + esc(r.note) : "");
}

return {esc, fmtE, FIDELITY, firstClause, shortLabel, renderTerms, groupRoutes,
        protonsShed, speciesCell, energyParts, routeEnergy, routeRung, routeCaveats,
        routeWhen, routeHead};
})();
