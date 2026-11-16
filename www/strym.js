// strym.js -- drag and drop for every file / folder field (the .dz boxes made by pathUI()).
// A browser reveals only the NAME and SIZE of a dropped file, never its location. So the page sends those to R, which looks for
// the file on this computer (see locate_dropped() in drives.R). Small files that cannot be found are uploaded as a copy through
// the field's hidden file input (the "strym-upload" message below). Very large files are never copied.
(function () {
  var stash = {};                                            // field id -> File objects of the last drop
  function zone(t) { return t && t.closest ? t.closest('.dz') : null; }
  function hasFiles(e) { return e.dataTransfer && [].slice.call(e.dataTransfer.types || []).indexOf('Files') >= 0; }

  ['dragenter', 'dragover'].forEach(function (ev) {
    document.addEventListener(ev, function (e) {
      if (!hasFiles(e)) return;
      e.preventDefault();                                    // otherwise the browser would open the dropped file
      var z = zone(e.target);
      document.querySelectorAll('.dz.over').forEach(function (x) { if (x !== z) x.classList.remove('over'); });
      if (z) { z.classList.add('over'); e.dataTransfer.dropEffect = 'copy'; } else e.dataTransfer.dropEffect = 'none';
    });
  });
  document.addEventListener('dragleave', function (e) {
    if (!e.relatedTarget) document.querySelectorAll('.dz.over').forEach(function (x) { x.classList.remove('over'); });
  });
  document.addEventListener('drop', function (e) {
    if (!hasFiles(e)) return;
    e.preventDefault();
    document.querySelectorAll('.dz.over').forEach(function (x) { x.classList.remove('over'); });
    var z = zone(e.target); if (!z) return;
    var ns = z.getAttribute('data-ns'), fl = e.dataTransfer.files, items = e.dataTransfer.items;
    var names = [], sizes = [], dirs = [], files = [];
    for (var i = 0; i < fl.length; i++) {
      var ent = null;
      try { ent = items && items[i] && items[i].webkitGetAsEntry ? items[i].webkitGetAsEntry() : null; } catch (err) { ent = null; }
      var isdir = !!(ent && ent.isDirectory);
      names.push(fl[i].name); sizes.push(isdir ? 0 : fl[i].size); dirs.push(isdir);
      if (!isdir) files.push(fl[i]);
    }
    if (!names.length) return;
    stash[ns] = files;
    Shiny.setInputValue(ns + '-dropped', { ns: ns, names: names, sizes: sizes, dirs: dirs, t: Date.now() }, { priority: 'event' });
  });

  window.addEventListener('load', function () {
    Shiny.addCustomMessageHandler('strym-upload', function (m) {
      var want = [].concat(m.names), fs = (stash[m.ns] || []).filter(function (f) { return want.indexOf(f.name) >= 0; });
      var inp = document.getElementById(m.ns + '-up');
      if (!fs.length || !inp) return;
      var dt = new DataTransfer(); fs.forEach(function (f) { dt.items.add(f); });
      inp.files = dt.files; inp.dispatchEvent(new Event('change', { bubbles: true }));
    });
  });
})();

// ---- the Browse window (see filebrowser.R) ------------------------------------------------------------------------------------------------
// Everything here is the browser's own work (no round trip to R): clicking a folder or a place tells R where to go, clicking a file selects it
// (Ctrl / Shift for several), the filter box hides rows as you type.
(function () {
  function send(id, v) { Shiny.setInputValue(id, v, { priority: 'event' }); }
  document.addEventListener('click', function (e) {
    var pl = e.target.closest && e.target.closest('.bf-place');
    if (pl) { send(pl.getAttribute('data-input'), { kind: pl.getAttribute('data-kind'), path: pl.getAttribute('data-path'), letter: pl.getAttribute('data-letter'), t: Date.now() }); return; }
    var list = e.target.closest && e.target.closest('.bf-list'); if (!list) return;
    var row = e.target.closest('.bf-row'); if (!row) return;
    if (row.classList.contains('bf-dir')) { send(list.getAttribute('data-cd'), list.getAttribute('data-base').replace(/\/+$/, '') + '/' + row.getAttribute('data-name')); return; }
    if (!row.classList.contains('bf-file')) return;
    var multi = list.getAttribute('data-multi') === 'true', rows = [].slice.call(list.querySelectorAll('.bf-file'));
    if (!multi) { rows.forEach(function (r) { r.classList.remove('sel'); }); row.classList.add('sel'); }
    else if (e.shiftKey && list._last) { var a = rows.indexOf(list._last), b = rows.indexOf(row); rows.slice(Math.min(a, b), Math.max(a, b) + 1).forEach(function (r) { r.classList.add('sel'); }); }
    else row.classList.toggle('sel');
    list._last = row;
    send(list.getAttribute('data-sel'), rows.filter(function (r) { return r.classList.contains('sel'); }).map(function (r) { return r.getAttribute('data-name'); }));
  });
  document.addEventListener('input', function (e) {
    var f = e.target; if (!f.classList || !f.classList.contains('bf-filter')) return;
    var q = f.value.toLowerCase();
    document.querySelectorAll('.bf-list .bf-row').forEach(function (r) { r.style.display = (!q || r.getAttribute('data-name').toLowerCase().indexOf(q) >= 0) ? '' : 'none'; });
  });
  document.addEventListener('keydown', function (e) {                         // Enter in the path box = Go
    if (e.key !== 'Enter' || !e.target.closest) return;
    var bar = e.target.closest('.bf-bar'); if (!bar) return;
    var b = bar.querySelector('button[id$="_go"]'); if (!b) return;
    e.preventDefault(); e.target.blur(); setTimeout(function () { b.click(); }, 350);
  });
})();
