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
