/* Keep the upstream engine, keyboard interaction and overlay. Benchmark metadata
 * adds descriptions and multiword matching without a parallel hand-written index. */
(() => {
  const base = new URL('../', document.currentScript.src);
  const original = Search.search;
  const render = Search.renderResults;
  const init = Search.init;
  const dispatch = Search.searchAndRender;
  let currentQuery = '', revision = 0;
  const canonical = value => value.toLowerCase().replace(/::/g, '/');
  Search.init = function(...args) {
    const result = init.apply(this, args);
    if (result !== false) {
      this.symbolCount = this.symbolCount.replace('symbols', 'benchmark pages');
      document.getElementById('search-symbolcount').textContent = this.symbolCount;
    }
    return result;
  };
  Search.search = function(query) {
    currentQuery = query.trim();
    const normal = canonical(query.trim());
    const [raw, completion] = original.call(this, query);
    if (!this.trie || !normal) return [raw, completion];
    const catalogue = window.RodinBenchmarkSearch || [];
    const terms = normal.split(/\s+/);
    const matches = catalogue.filter(item => {
      const text = canonical(item.path + ' ' + item.display + ' ' + item.title + ' ' + item.target + ' ' + item.text + ' ' + item.brief);
      return terms.every(term => text.includes(term));
    }).sort((a, b) => {
      const score = item => canonical(item.title) === normal ? 0 : canonical(item.title).startsWith(normal) ? 1 : canonical(item.path).includes(normal) ? 2 : 3;
      return score(a) - score(b) || a.path.split('/').length - b.path.split('/').length || a.display.localeCompare(b.display);
    });
    const saved = { string: this.searchString, stack: this.searchStack.slice() };
    const results = [];
    const seen = new Set();
    try {
      for (const entry of matches) {
        const [found] = original.call(this, entry.path);
        const result = found.find(item => item.url === entry.url);
        if (!result || seen.has(entry.url)) continue;
        seen.add(entry.url);
        const name = this.fromUtf8(result.name);
        const position = name.toLowerCase().indexOf(query.trim().toLowerCase());
        results.push({ ...result, url: new URL(entry.url, base).href,
          suffixLength: position < 0 ? 0 : this.toUtf8(name.slice(position + query.trim().length)).length });
        if (results.length === 50) break;
      }
    } finally {
      this.searchString = saved.string;
      this.searchStack = saved.stack;
    }
    // Only the unchanged prefix result set can offer meaningful completion.
    return [results, terms.length === 1 && results.length && raw.length ? completion : ''];
  };
  Search.searchAndRender = function(value) {
    dispatch.call(this, value);
    if (value.trim() && !this.searchString.length)
      document.getElementById('search-symbolcount').textContent = document.querySelectorAll('#search-results > li').length + ' results';
  };
  Search.renderResults = function(results) {
    const current = ++revision;
    const list = document.getElementById('search-results');
    window.MathJax?.typesetClear?.([list]);
    // The engine retains only the prefix present in its trie. Descriptions and
    // multiword queries still need the full user input for overlay state.
    const saved = this.searchString;
    try {
      this.searchString = this.toUtf8(currentQuery.toLowerCase());
      render.call(this, results);
    } finally { this.searchString = saved; }
    const catalogue = window.RodinBenchmarkSearch || [];
    for (const link of document.querySelectorAll('#search-results > li > a')) {
      const item = catalogue.find(item => new URL(item.url, base).href === link.href);
      if (!item) continue;
      const name = link.lastElementChild;
      name.className = 'benchmark-search-path';
      name.textContent = item.display;
      link.dataset.mdLinkTitle = item.title;
      const brief = document.createElement('p');
      brief.className = 'benchmark-search-brief';
      brief.textContent = item.brief;
      link.append(brief);
      const terms = canonical(currentQuery).split(/[\s/]+/).filter(Boolean);
      for (const node of [name, brief]) {
        const text = node.textContent;
        node.replaceChildren();
        // Keep TeX delimiters and commands intact for the shared math renderer.
        for (const part of text.split(/(\\\([\s\S]*?\\\)|\\\[[\s\S]*?\\\])/g)) {
          if (part.startsWith('\\(') || part.startsWith('\\[')) { node.append(document.createTextNode(part));continue; }
          let start = 0;
          while (start < part.length) {
            const matches = terms.map(term => ({ term, at: part.toLowerCase().indexOf(term, start) })).filter(match => match.at >= 0).sort((a,b) => a.at - b.at || b.term.length - a.term.length);
            if (!matches.length) { node.append(document.createTextNode(part.slice(start)));break; }
            const match = matches[0];
            node.append(document.createTextNode(part.slice(start, match.at)));
            const typed = document.createElement('span');typed.className = 'm-doc-search-typed';typed.textContent = part.slice(match.at, match.at + match.term.length);node.append(typed);
            start = match.at + match.term.length;
          }
        }
      }
    }
    if (window.MathJax?.startup?.promise) {
      this.mathReady = window.MathJax.startup.promise.then(() => {
        if (current === revision) return window.MathJax.typesetPromise([list]);
      });
    }
  };
})();

// Keep the documentation shortcuts from consuming text entered in chart filters
// or the theme menu. Closing the overlay restores the invoking control's focus.
(() => {
  const keyboard = document.onkeydown;
  document.onkeydown = function(event) {
    if (location.hash !== '#search' && (event.defaultPrevented || event.target.closest?.('input, textarea, select, [contenteditable="true"], [role="menu"]'))) return;
    return keyboard.call(this, event);
  };
  let previous;
  const open = window.showSearch, close = window.hideSearch;
  window.showSearch = function() { previous = document.activeElement;return open(); };
  window.hideSearch = function() {
    const result = close();
    if (previous?.isConnected && !previous.closest('#search')) previous.focus();
    else document.querySelector('.benchmark-search-trigger')?.focus();
    return result;
  };
})();
