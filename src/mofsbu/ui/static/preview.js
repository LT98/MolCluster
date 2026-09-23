"use strict";
/* Showing a structure's geometry without leaving the page.

   Two shapes, because there are two questions.  A popover answers "which one is that?"
   while the pointer is passing over a name, and disappears.  An embedded viewer answers
   "what does this actually look like?" about the one structure the page is currently
   about, and stays, and can be turned around.

   **One viewer each, never one per use.** A browser caps live WebGL contexts at eight to
   sixteen; a context per hovered link would reach that in a dozen names and then start
   losing the earliest ones. So the popover keeps a single viewer and swaps the model in
   it, and `embed` hands back one viewer per host element.

   Coordinates come from `/api/structures/{id}/xyz`, which exists precisely because the
   full structure endpoint prices every incoming route and is far too dear at hover rate.
   A structure with no stored geometry is a fact about the registry, so it is shown as
   one rather than left as an empty box. */

const mofsbuPreview = (() => {

const CACHE = new Map();
const dark = () => matchMedia("(prefers-color-scheme: dark)").matches;
const bg = () => dark() ? "#16181c" : "#f7f7f8";
const has3D = () => typeof $3Dmol !== "undefined";

async function xyzOf(sid){
  if(CACHE.has(sid)) return CACHE.get(sid);
  const r = await fetch(`/api/structures/${sid}/xyz`);
  if(!r.ok) throw new Error(r.status === 404
    ? "no geometry stored for this structure" : `could not load structure ${sid}`);
  const text = await r.text();
  CACHE.set(sid, text);
  return text;
}

function paint(viewer, xyz, style){
  viewer.removeAllModels();
  viewer.addModel(xyz, "xyz");
  viewer.setStyle({}, STYLES[style] || STYLES.stick);
  viewer.zoomTo();
  viewer.render();
}

const STYLES = {
  stick: {stick: {radius: 0.13}, sphere: {scale: 0.18}},
  ball:  {stick: {radius: 0.10}, sphere: {scale: 0.30}},
  space: {sphere: {scale: 0.90}},
};

// ── the popover ────────────────────────────────────────────────────────────
let pop = null, popViewer = null, timer = null, shown = null;

function mount(){
  if(pop) return pop;
  pop = document.createElement("div");
  pop.id = "mofsbu-preview";
  pop.hidden = true;
  pop.innerHTML = `<div class="mp-3d"></div><div class="mp-note"></div>`;
  document.body.appendChild(pop);
  return pop;
}

function place(el){
  const r = el.getBoundingClientRect(), box = pop.getBoundingClientRect();
  const gap = 10;
  let left = r.right + gap;
  if(left + box.width > innerWidth - 6) left = Math.max(6, r.left - box.width - gap);
  let top = r.top + r.height / 2 - box.height / 2;
  top = Math.max(6, Math.min(top, innerHeight - box.height - 6));
  pop.style.left = `${Math.round(left)}px`;
  pop.style.top = `${Math.round(top)}px`;
}

async function show(el, sid, caption){
  mount();
  shown = sid;
  pop.hidden = false;
  pop.querySelector(".mp-note").textContent = caption || `structure ${sid}`;
  place(el);
  if(!has3D()){
    pop.querySelector(".mp-3d").textContent = "3Dmol did not load";
    return;
  }
  let xyz;
  try { xyz = await xyzOf(sid); }
  catch(e){
    if(shown !== sid) return;
    pop.querySelector(".mp-3d").textContent = e.message;
    place(el);
    return;
  }
  if(shown !== sid) return;                 // the pointer moved on while we fetched
  const host = pop.querySelector(".mp-3d");
  if(!popViewer){
    host.textContent = "";
    popViewer = $3Dmol.createViewer(host, {backgroundColor: bg()});
  }
  paint(popViewer, xyz, "stick");
  place(el);
}

function hide(){
  clearTimeout(timer);
  shown = null;
  if(pop) pop.hidden = true;
}

/** Hover any `[data-sid]` inside `container` for `delay` ms and its geometry appears. */
function bind(container, opts){
  const {delay = 450, caption} = opts || {};
  container.addEventListener("mouseover", (e) => {
    const el = e.target.closest("[data-sid]");
    if(!el || !el.dataset.sid) return;
    clearTimeout(timer);
    const sid = Number(el.dataset.sid);
    timer = setTimeout(() => show(el, sid, caption && caption(sid, el)), delay);
  });
  container.addEventListener("mouseout", (e) => {
    if(e.target.closest("[data-sid]")) hide();
  });
  // A popover pinned to an element that has scrolled away is worse than none.
  container.addEventListener("scroll", hide, {passive: true});
}

// ── the embedded viewer ────────────────────────────────────────────────────
/** One persistent, turnable viewer inside `host`. Returns `{show, clear, style}`. */
function embed(host){
  let viewer = null, current = null, style = "stick";
  const note = document.createElement("div");
  note.className = "mv-note";
  host.appendChild(note);
  const stage = document.createElement("div");
  stage.className = "mv-stage";
  host.appendChild(stage);

  function say(text){ note.textContent = text || ""; note.hidden = !text; }

  return {
    async show(sid){
      if(sid === current) return;
      current = sid;
      if(sid == null){ say("nothing selected"); stage.hidden = true; return; }
      if(!has3D()){ say("3Dmol did not load"); stage.hidden = true; return; }
      let xyz;
      try { xyz = await xyzOf(sid); }
      catch(e){ if(current === sid){ say(e.message); stage.hidden = true; } return; }
      if(current !== sid) return;
      stage.hidden = false;
      say("");
      if(!viewer) viewer = $3Dmol.createViewer(stage, {backgroundColor: bg()});
      paint(viewer, xyz, style);
      viewer.resize();
    },
    style(next){
      style = next;
      if(viewer){ viewer.setStyle({}, STYLES[next] || STYLES.stick); viewer.render(); }
    },
    recentre(){ if(viewer){ viewer.zoomTo(); viewer.render(); } },
    resize(){ if(viewer) viewer.resize(); },
  };
}

return {mount, bind, embed, hide, xyzOf, STYLES};
})();
