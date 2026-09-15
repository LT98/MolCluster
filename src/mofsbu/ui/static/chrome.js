"use strict";
/* Chrome the three pages share: the tab strip, and the control that says which
   registry the server is pointed at.

   Both are here rather than in each page because both are properties of the SERVER,
   not of a page.  There is one list of pages, so a new one appears in all three tab
   strips at once; and there is one database switch, so the viewer and the builder
   cannot disagree about which registry is current (ui/active.py explains why that
   matters). */

/* The tab strip, in order.  Adding a page is this line and its route — nothing in
   any of the three HTML files changes.  `soon:true` reserves the space for a page
   that is named but not built; it renders dimmed and does not navigate. */
const MOFSBU_PAGES = [
  {href: "/",        label: "registry",       title: "browse stored structures"},
  {href: "/builder", label: "builder",        title: "write a spec and queue a run"},
  {href: "/runs",    label: "run inspector",  title: "what each run planned, built and refused"},
];

const mofsbuChrome = {
  esc(s){ return String(s ?? "").replace(/[&<>"]/g, c =>
    ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c])); },

  /** A page's identity for comparison: "/runs/" and "/runs" are the same page. */
  norm(href){ return href.replace(/\/+$/, "") || "/"; },

  mountTabs(host){
    const now = this.norm(location.pathname);
    if(!host.getAttribute("aria-label")) host.setAttribute("aria-label", "sections");
    host.innerHTML = MOFSBU_PAGES.map(p => {
      const current = this.norm(p.href) === now;
      if(p.soon)
        return `<span class="tab soon" title="${this.esc(p.title || "not built yet")}"`
             + ` aria-disabled="true">${this.esc(p.label)}</span>`;
      return `<a class="tab" href="${this.esc(p.href)}" title="${this.esc(p.title || "")}"`
           + `${current ? ' aria-current="page"' : ""}>${this.esc(p.label)}</a>`;
    }).join("");
  },
};

/* ── the active registry ──────────────────────────────────────────────────────
   Wraps /api/databases and /api/database/select so the viewer and the builder ask
   the same questions the same way.  The awkward part is labelling: two folders can
   hold a `registry.db` and the launch database may sit outside the data root, so the
   option VALUE is always a full path and the label disambiguates only when it has
   to.  A dropdown with two identical entries is not a choice. */
const mofsbuDb = {
  async _json(url, opts){
    const r = await fetch(url, opts);
    // These endpoints refuse with a sentence in `detail` — 404 on a path the server
    // never offered, 409 on an ambiguous name.  Showing that beats showing the JSON.
    const body = await r.json().catch(() => ({}));
    if(!r.ok) throw new Error(body.detail || body.error || r.statusText);
    return body;
  },
  list(){ return this._json("/api/databases"); },
  select(path){
    return this._json("/api/database/select", {
      method: "POST", headers: {"Content-Type": "application/json"},
      body: JSON.stringify({path})});
  },
  create(name){
    return this._json("/api/databases", {
      method: "POST", headers: {"Content-Type": "application/json"},
      body: JSON.stringify({name})});
  },

  /** The row count is the point: "which one is my real one" is the question being
      asked, and two plausible filenames do not answer it. */
  describe(d){
    return d.error ? "unreadable"
         : !d.exists ? "empty (not created yet)"
         : `${d.n_structures} structures · ${d.n_runs} runs`;
  },

  fill(select, body){
    const esc = mofsbuChrome.esc;
    const counts = {};
    body.databases.forEach(d => counts[d.name] = (counts[d.name] || 0) + 1);
    select.innerHTML = body.databases.map(d => {
      const where = counts[d.name] > 1 ? ` [${d.folder}]` : "";
      return `<option value="${esc(d.path)}" title="${esc(d.path)}"`
           + `${d.active ? " selected" : ""}>`
           + `${esc(d.name + where)} — ${esc(this.describe(d))}</option>`;
    }).join("");
    return body.databases.find(d => d.active) || {};
  },

  /** Populate a <select>, and switch the whole app when it changes.
      `noteText(active, body)` is the line to show under it; `onSwitch(path)` is what
      the calling page has to redo once the registry underneath it has moved. */
  async mount(opts){
    const {select, note, onSwitch, onError, noteText} = opts;
    let body;
    try { body = await this.list(); }
    catch(e){ if(note) note.textContent = `could not list databases — ${e.message}`; return null; }
    const active = this.fill(select, body);
    if(note && noteText) note.textContent = noteText(active, body);
    select.onchange = async () => {
      const to = select.value;
      try { await this.select(to); }
      catch(e){
        // The server kept the registry it had, so re-list rather than leave the
        // dropdown showing a switch that did not happen.
        if(onError) onError(e); else if(note) note.textContent = e.message;
        await this.mount(opts);
        return;
      }
      await this.mount(opts);
      if(onSwitch) await onSwitch(to);
    };
    return body;
  },
};

document.addEventListener("DOMContentLoaded", () => {
  document.querySelectorAll("nav.tabs").forEach(n => mofsbuChrome.mountTabs(n));
});
