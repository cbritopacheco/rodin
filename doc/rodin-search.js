/* Rank general class references before specializations for class-name queries.
 * The catalogue is filled from Doxygen XML by documentation_search.py.
 */
(function(search, catalogue) {
  'use strict';

  const original = search.search;
  search.search = function(query) {
    const normalized = query.trim().toLowerCase();
    // Template arguments and scope-list requests express a more specific intent.
    if(!this.trie || !normalized || normalized.includes('<') || /(?:::|\/)$/.test(normalized))
      return original.call(this, query);

    const matches = catalogue.map(entry => {
      const key = entry.keys.find(key => key.startsWith(normalized));
      return key ? {entry, remainder: key.length - normalized.length} : null;
    }).filter(Boolean).sort((a, b) =>
      a.remainder - b.remainder || a.entry.name.length - b.entry.name.length ||
      a.entry.name.localeCompare(b.entry.name));
    if(!matches.length)
      return original.call(this, query);

    const savedString = this.searchString;
    const savedStack = this.searchStack.slice();
    const promoted = [];
    const seen = new Set();
    try {
      for(const match of matches) {
        // Query the existing index rather than constructing synthetic results.
        // This also finds primaries excluded by the broad query's result limit.
        const [results] = original.call(this, match.entry.name);
        const result = results.find(result => result.url === match.entry.url);
        if(result && !(result.flags & 6) && !seen.has(result.url)) {
          promoted.push(result);
          seen.add(result.url);
          if(promoted.length >= this.maxResults)
            break;
        }
      }
    } finally {
      // m.css caches the last query's trie path. Restore the user's query even
      // after a failed lookup, along with its completion and highlighting state.
      this.searchString = savedString;
      this.searchStack = savedStack;
    }
    const [results, completion] = original.call(this, query);
    // Exact lookups have different highlight suffixes. A returned primary must
    // highlight the user's prefix, not pretend the complete name was typed.
    for(const result of promoted) {
      const existing = results.find(candidate => candidate.url === result.url);
      if(existing) {
        result.suffixLength = existing.suffixLength;
      } else {
        const name = this.toUtf8(this.fromUtf8(result.name).toLowerCase());
        // Result names are UTF-8 byte strings in the upstream search API.
        const encoded = this.toUtf8(normalized);
        const index = name.lastIndexOf(encoded);
        result.suffixLength = index < 0 ? 0 : name.length - index - encoded.length;
      }
    }
    return [promoted.concat(results.filter(result => !seen.has(result.url))).slice(0, this.maxResults), completion];
  };
})(Search, __RODIN_SEARCH_CLASSES__);
