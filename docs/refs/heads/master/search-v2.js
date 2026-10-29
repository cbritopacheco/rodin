/*
    This file is part of m.css.

    Copyright © 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025
              Vladimír Vondruš <mosra@centrum.cz>
    Copyright © 2020 Sergei Izmailov <sergei.a.izmailov@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included
    in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.
*/

"use strict"; /* it summons the Cthulhu in a proper way, they say */

var Search = {
    formatVersion: 2, /* the data filename contains this number too */

    dataSize: 0, /* used mainly by tests, not here */
    symbolCount: '&hellip;',
    trie: null,
    map: null,
    mapFlagsOffset: null,
    typeMap: null,
    maxResults: 0,

    /* Type sizes and masks. The data is always fetched as 16/32bit number and
       then masked to 1, 2, 3 or 4 bytes. Fortunately on LE a mask is enough,
       on BE we'd have to read N bytes before and then mask. */
    nameSizeBytes: null,
    nameSizeMask: null,
    resultIdBytes: null,
    resultIdMask: null,
    fileOffsetBytes: null,
    fileOffsetMask: null,
    lookaheadBarrierMask: null,

    /* Always contains at least the root node offset and then one node offset
       per entered character */
    searchString: '',
    searchStack: [],

    /* So items don't get selected right away when a cursor is over results but
       only after mouse moves */
    mouseMovedSinceLastRender: false,

    /* Whether we can go back in history in order to hide the search box or
       not. We can't do that if we arrived directly on #search from outside. */
    canGoBackToHideSearch: false,

    /* Autocompletion in the input field is whitelisted only for character
       input (so not deletion, cut, or anything else). This is flipped in the
       onkeypress event and reset after each oninput event. */
    autocompleteNextInputEvent: false,

    init: function(buffer, maxResults) {
        let view = new DataView(buffer);

        /* The file is too short to contain at least the headers and empty
           sections */
        if(view.byteLength < 31) {
            console.error("Search data too short");
            return false;
        }

        if(view.getUint8(0) != 'M'.charCodeAt(0) ||
           view.getUint8(1) != 'C'.charCodeAt(0) ||
           view.getUint8(2) != 'S'.charCodeAt(0)) {
            console.error("Invalid search data signature");
            return false;
        }

        if(view.getUint8(3) != this.formatVersion) {
            console.error("Invalid search data version");
            return false;
        }

        /* Fetch type sizes. The only value that can fail is result ID byte
           count, where value of 3 has no assigned meaning. */
        let typeSizes = view.getUint8(4, true);
        if((typeSizes & 0x01) >> 0 == 0) {
            this.fileOffsetBytes = 3;
            this.fileOffsetMask = 0x00ffffff;
            this.lookaheadBarrierMask = 0x00800000;
        } else /* (typeSizes & 0x01) >> 0 == 1 */ {
            this.fileOffsetBytes = 4;
            this.fileOffsetMask = 0xffffffff;
            this.lookaheadBarrierMask = 0x80000000;
        }
        if((typeSizes & 0x06) >> 1 == 0) {
            this.resultIdBytes = 2;
            this.resultIdMask = 0x0000ffff;
        } else if((typeSizes & 0x06) >> 1 == 1) {
            this.resultIdBytes = 3;
            this.resultIdMask = 0x00ffffff;
        } else if((typeSizes & 0x06) >> 1 == 2) {
            this.resultIdBytes = 4;
            this.resultIdMask = 0xffffffff;
        } else /* (typeSizes & 0x06) >> 1 == 3 */ {
            console.error("Invalid search data result ID byte value");
            return false;
        }
        if((typeSizes & 0x08) >> 3 == 0) {
            this.nameSizeBytes = 1;
            this.nameSizeMask = 0x00ff;
        } else /* (typeSizes & 0x08) >> 3 == 1 */ {
            this.nameSizeBytes = 2;
            this.nameSizeMask = 0xffff;
        }

        /* Separate the data into the trie and the result / type map. Because
           we're reading larger values than there might be and then masking out
           the high bytes, keep extra 1/2 byte padding at the end to avoid
           OOB errors. */
        let mapOffset = view.getUint32(12, true);
        let typeMapOffset = view.getUint32(16, true);
        /* There may be a 3-byte file offset at the end of the trie which we'll
           read as 32-bit, add one safety byte in that case */
        this.trie = new DataView(buffer, 20, mapOffset - 20 + (4 - this.fileOffsetBytes));
        /* There may be a 3-byte file size (for zero results) which we'll read
           as 32-bit, add one safety byte in that case */
        this.map = new DataView(buffer, mapOffset, typeMapOffset - mapOffset + (4 - this.fileOffsetBytes));
        /* No variable-size types in the type map at the moment */
        this.typeMap = new DataView(buffer, typeMapOffset);

        /* Offset of the first result map item is after N + 1 offsets and N
           flags, calculate flag offset from that */
        this.mapFlagsOffset = this.fileOffsetBytes*(((this.map.getUint32(0, true) & this.fileOffsetMask) - this.fileOffsetBytes)/(this.fileOffsetBytes + 1) + 1);

        /* Set initial properties */
        this.dataSize = buffer.byteLength;
        this.symbolCount = view.getUint32(8, true) + " symbols (" + Math.round(this.dataSize/102.4)/10 + " kB)";
        this.maxResults = maxResults ? maxResults : 100;
        this.searchString = '';
        this.searchStack = [this.trie.getUint32(0, true)];

        /* istanbul ignore if */
        if(typeof document !== 'undefined') {
            document.getElementById('search-symbolcount').innerHTML = this.symbolCount;
            document.getElementById('search-input').disabled = false;
            document.getElementById('search-input').placeholder = "Type something here …";
            document.getElementById('search-input').focus();

            /* Search for the input value (there might be something already,
               for example when going back in the browser) */
            let value = document.getElementById('search-input').value;

            /* Otherwise check the GET parameters for `q` and fill the input
               with that */
            if(!value.length) {
                var args = decodeURIComponent(window.location.search.substr(1)).trim().split('&');
                for(var i = 0; i != args.length; ++i) {
                    if(args[i].substring(0, 2) != 'q=') continue;

                    value = document.getElementById('search-input').value = args[i].substring(2);
                    break;
                }
            }

            if(value.length) Search.searchAndRender(value);
        }

        return true;
    },

    download: /* istanbul ignore next */ function(url) {
        var req = window.XDomainRequest ? new XDomainRequest() : new XMLHttpRequest();
        if(!req) return;

        req.open("GET", url, true);
        req.responseType = 'arraybuffer';
        req.onreadystatechange = function() {
            if(req.readyState != 4) return;

            Search.init(req.response);
        }
        req.send();
    },

    base85decode: function(base85string) {
        function charValue(char) {
            if(char >=  48 && char <  58) /* 0-9 -> 0-9 */
                return char - 48 + 0;
            if(char >=  65 && char <  91) /* A-Z -> 10-35 */
                return char - 65 + 10;
            if(char >=  97 && char < 123) /* a-z -> 36-61 */
                return char - 97 + 36;
            if(char ==  33)               /*  !  -> 62 */
                return 62;
            /* skipping 34 (') */
            if(char >=  35 && char <  39) /* #-& -> 63-66 */
                return char - 35 + 63;
            /* skipping 39 (") */
            if(char >=  40 && char <  44) /* (-+ -> 67-70 */
                return char - 40 + 67;
            /* skipping 44 (,) */
            if(char ==  45)               /*  -  -> 71 */
                return 71;
            if(char >=  59 && char <  65) /* ;-@ -> 72-77 */
                return char - 59 + 72;
            if(char >=  94 && char <  97) /* ^-` -> 78-80 */
                return char - 94 + 78;
            if(char >= 123 && char < 127) /* {-~ -> 81-84 */
                return char - 123 + 81;

            return 0; /* Interpret padding values as zeros */
        }

        /* Pad the string for easier decode later. We don't read past the file
           end, so it doesn't matter what garbage is there. */
        if(base85string.length % 5) {
            console.log("Expected properly padded base85 data");
            return;
        }

        let buffer = new ArrayBuffer(base85string.length*4/5);
        let data8 = new DataView(buffer);
        for(let i = 0; i < base85string.length; i += 5) {
            let char1 = charValue(base85string.charCodeAt(i + 0));
            let char2 = charValue(base85string.charCodeAt(i + 1));
            let char3 = charValue(base85string.charCodeAt(i + 2));
            let char4 = charValue(base85string.charCodeAt(i + 3));
            let char5 = charValue(base85string.charCodeAt(i + 4));

            data8.setUint32(i*4/5, char5 +
                                   char4*85 +
                                   char3*85*85 +
                                   char2*85*85*85 +
                                   char1*85*85*85*85, false); /* BE, yes */
        }

        return buffer;
    },

    load: function(base85string) {
        return this.init(this.base85decode(base85string));
    },

    /* http://ecmanaut.blogspot.com/2006/07/encoding-decoding-utf8-in-javascript.html */
    toUtf8: function(string) { return unescape(encodeURIComponent(string)); },
    fromUtf8: function(string) { return decodeURIComponent(escape(string)); },

    autocompletedCharsToUtf8: function(chars) {
        /* Strip incomplete UTF-8 chars from the autocompletion end */
        for(let i = chars.length - 1; i >= 0; --i) {
            let c = chars[i];

            /* We're safe, finish */
            if(
                /* ASCII value at the end */
                (c < 128 && i + 1 == chars.length) ||

                /* Full two-byte character at the end */
                ((c & 0xe0) == 0xc0 && i + 2 == chars.length) ||

                /* Full three-byte character at the end */
                ((c & 0xf0) == 0xe0 && i + 3 == chars.length) ||

                /* Full four-byte character at the end */
                ((c & 0xf8) == 0xf0 && i + 4 == chars.length)
            ) break;

            /* Continuing UTF-8 character, go further back */
            if((c & 0xc0) == 0x80) continue;

            /* Otherwise the character is not complete, drop it from the end */
            chars.length = i;
            break;
        }

        /* Convert the autocompleted UTF-8 sequence to a string */
        let suggestedTabAutocompletionString = '';
        for(let i = 0; i != chars.length; ++i)
            suggestedTabAutocompletionString += String.fromCharCode(chars[i]);
        return suggestedTabAutocompletionString;
    },

    /* Returns the values in UTF-8, but input is in whatever shitty 16bit
       encoding JS has */
    search: function(searchString) {
        /* Normalize the search string first, convert to UTF-8 and trim spaces
           from the left. From the right they're trimmed only if nothing is
           found, see below. */
        searchString = this.toUtf8(searchString.toLowerCase().replace(/^\s+/,''));

        /* TODO: maybe i could make use of InputEvent.data and others here */

        /* Find longest common prefix of previous and current value so we don't
           need to needlessly search again */
        let max = Math.min(searchString.length, this.searchString.length);
        let commonPrefix = 0;
        for(; commonPrefix != max; ++commonPrefix)
            if(searchString[commonPrefix] != this.searchString[commonPrefix]) break;

        /* Drop items off the stack if it has has more than is needed for the
           common prefix (it needs to have at least one item, though) */
        if(commonPrefix + 1 < this.searchStack.length)
            this.searchStack.splice(commonPrefix + 1, this.searchStack.length - commonPrefix - 1);

        /* Add new characters from the search string */
        let foundPrefix = commonPrefix;
        for(; foundPrefix != searchString.length; ++foundPrefix) {
            /* Calculate offset and count of children */
            let offset = this.searchStack[this.searchStack.length - 1];

            /* If there's a lot of results, the result count is a 16bit BE value
               instead */
            let resultCount = this.trie.getUint8(offset);
            let resultCountSize = 1;
            if(resultCount & 0x80) {
                resultCount = this.trie.getUint16(offset, false) & ~0x8000;
                ++resultCountSize;
            }

            let childCount = this.trie.getUint8(offset + resultCountSize);

            /* Go through all children and find the next offset */
            let childOffset = offset + resultCountSize + 1 + resultCount*this.resultIdBytes;
            let found = false;
            for(let j = 0; j != childCount; ++j) {
                if(String.fromCharCode(this.trie.getUint8(childOffset + j)) != searchString[foundPrefix])
                    continue;

                this.searchStack.push(this.trie.getUint32(childOffset + childCount + j*this.fileOffsetBytes, true) & this.fileOffsetMask & ~this.lookaheadBarrierMask);
                found = true;
                break;
            }

            /* Character not found */
            if(!found) {
                /* If we found everything except spaces at the end, pretend the
                   spaces aren't there. On the other hand, we *do* want to
                   try searching with the spaces first -- it can narrow down
                   the result list for page names or show subpages (which are
                   after a lookahead barrier that's a space). */
                if(!searchString.substr(foundPrefix).trim().length)
                    searchString = searchString.substr(0, foundPrefix);

                break;
            }
        }

        /* Save the whole found prefix for next time */
        this.searchString = searchString.substr(0, foundPrefix);

        /* If the whole thing was not found, return an empty result and offer
           external search */
        if(foundPrefix != searchString.length) {
            /* istanbul ignore if */
            if(typeof document !== 'undefined') {
                let link = document.getElementById('search-external');
                if(link)
                    link.href = link.dataset.searchEngine.replace('{query}', encodeURIComponent(searchString));
            }
            return [[], ''];
        }

        /* Otherwise gather the results */
        let suggestedTabAutocompletionChars = [];
        let results = [];
        let leaves = [[this.searchStack[this.searchStack.length - 1], 0]];
        while(leaves.length) {
            /* Pop offset from the queue */
            let current = leaves.shift();
            let offset = current[0];
            let suffixLength = current[1];

            /* Calculate child count. If there's a lot of results, the count
               "leaks over" to the child count storage. */
            /* TODO: hmmm. this is helluvalot duplicated code. hmm. */
            let resultCount = this.trie.getUint8(offset);
            let resultCountSize = 1;
            if(resultCount & 0x80) {
                resultCount = this.trie.getUint16(offset, false) & ~0x8000;
                ++resultCountSize;
            }

            let childCount = this.trie.getUint8(offset + resultCountSize);

            /* Populate the results with all values associated with this node */
            for(let i = 0; i != resultCount; ++i) {
                let index = this.trie.getUint32(offset + resultCountSize + 1 + i*this.resultIdBytes, true) & this.resultIdMask;
                results.push(this.gatherResult(index, suffixLength, 0xffffff)); /* should be enough haha */

                /* 'nuff said. */
                if(results.length >= this.maxResults)
                    return [results, this.autocompletedCharsToUtf8(suggestedTabAutocompletionChars)];
            }

            /* Dig deeper */
            let childOffset = offset + resultCountSize + 1 + resultCount*this.resultIdBytes;
            for(let j = 0; j != childCount; ++j) {
                let offsetBarrier = this.trie.getUint32(childOffset + childCount + j*this.fileOffsetBytes, true) & this.fileOffsetMask;

                /* Lookahead barrier, don't dig deeper */
                if(offsetBarrier & this.lookaheadBarrierMask) continue;

                /* Append to the queue */
                leaves.push([offsetBarrier & ~this.lookaheadBarrierMask, suffixLength + 1]);

                /* We don't have anything yet and this is the only path
                   forward, add the char to suggested Tab autocompletion. Can't
                   extract it from the leftmost 8 bits of offsetBarrier because
                   that would make it negative, have to load as Uint8 instead.
                   Also can't use String.fromCharCode(), because later doing
                   str.charCodeAt() would give me back UTF-16 values, which is
                   absolutely unwanted when all I want is check for truncated
                   UTF-8. */
                if(!results.length && leaves.length == 1 && childCount == 1)
                    suggestedTabAutocompletionChars.push(this.trie.getUint8(childOffset + j));
            }
        }

        return [results, this.autocompletedCharsToUtf8(suggestedTabAutocompletionChars)];
    },

    gatherResult: function(index, suffixLength, maxUrlPrefix) {
        let flags = this.map.getUint8(this.mapFlagsOffset + index);
        let resultOffset = this.map.getUint32(index*this.fileOffsetBytes, true) & this.fileOffsetMask;

        /* The result is an alias, parse the aliased prefix */
        let aliasedIndex = null;
        if((flags & 0xf0) == 0x00) {
            aliasedIndex = this.map.getUint32(resultOffset, true) & this.resultIdMask;
            resultOffset += this.resultIdBytes;
        }

        /* The result has a prefix, parse that first, recursively */
        let name = '';
        let url = '';
        if(flags & (1 << 3)) {
            let prefixIndex = this.map.getUint32(resultOffset, true) & this.resultIdMask;
            let prefixUrlPrefixLength = Math.min(this.map.getUint16(resultOffset + this.resultIdBytes, true) & this.nameSizeMask, maxUrlPrefix);

            let prefix = this.gatherResult(prefixIndex, 0 /*ignored*/, prefixUrlPrefixLength);
            name = prefix.name;
            url = prefix.url;

            resultOffset += this.resultIdBytes + this.nameSizeBytes;
        }

        /* The result has a suffix, extract its length */
        let resultSuffixLength = 0;
        if(flags & (1 << 0)) {
            resultSuffixLength = this.map.getUint16(resultOffset, true) & this.nameSizeMask;
            resultOffset += this.nameSizeBytes;
        }

        let nextResultOffset = this.map.getUint32((index + 1)*this.fileOffsetBytes, true) & this.fileOffsetMask;

        /* Extract name */
        let j = resultOffset;
        for(; j != nextResultOffset; ++j) {
            let c = this.map.getUint8(j);

            /* End of null-delimited name */
            if(!c) {
                ++j;
                break; /* null-delimited */
            }

            name += String.fromCharCode(c); /* eheh. IS THIS FAST?! */
        }

        /* The result is an alias and we're not deep inside resolving a prefix,
           extract the aliased name and URL */
        /* TODO: this abuses 0xffffff to guess how the call stack is deep and
           that's just wrong, fix! */
        if(aliasedIndex != null && maxUrlPrefix == 0xffffff) {
            let alias = this.gatherResult(aliasedIndex, 0 /* ignored */, 0xffffff); /* should be enough haha */

            /* Keeping in UTF-8, as we need that for proper slicing (and concatenating) */
            return {name: name,
                    alias: alias.name,
                    url: alias.url,
                    flags: alias.flags,
                    cssClass: alias.cssClass,
                    typeName: alias.typeName,
                    suffixLength: suffixLength + resultSuffixLength};
        }

        /* Otherwise extract URL from here */
        let max = Math.min(j + maxUrlPrefix - url.length, nextResultOffset);
        for(; j != max; ++j) {
            url += String.fromCharCode(this.map.getUint8(j));
        }

        /* This is an alias, return what we have, without parsed CSS class and
           type name as those are retrieved from the final target type */
        if(!(flags >> 4))
            return {name: name,
                    url: url,
                    flags: flags & 0x0f,
                    suffixLength: suffixLength + resultSuffixLength};

        /* Otherwise, get CSS class and type name for the result label */
        let typeMapIndex = (flags >> 4) - 1;
        let cssClass = [
            /* Keep in sync with _search.py */
            'm-default',
            'm-primary',
            'm-success',
            'm-warning',
            'm-danger',
            'm-info',
            'm-dim'
        ][this.typeMap.getUint8(typeMapIndex*2)];
        let typeNameOffset = this.typeMap.getUint8(typeMapIndex*2 + 1);
        let nextTypeNameOffset = this.typeMap.getUint8((typeMapIndex + 1)*2 + 1);
        let typeName = '';
        for(let j = typeNameOffset; j != nextTypeNameOffset; ++j)
            typeName += String.fromCharCode(this.typeMap.getUint8(j));

        /* Keeping in UTF-8, as we need that for proper slicing (and
           concatenating). Strip the type from the flags, as it's now expressed
           directly. */
        return {name: name,
                url: url,
                flags: flags & 0x0f,
                cssClass: cssClass,
                typeName: typeName,
                suffixLength: suffixLength + resultSuffixLength};
    },

    escape: function(name) {
        return name.replace(/[\"&<>]/g, function (a) {
            return { '"': '&quot;', '&': '&amp;', '<': '&lt;', '>': '&gt;' }[a];
        });
    },
    escapeForRtl: function(name) {
        /* Besides the obvious escaping of HTML entities we also need
           to escape punctuation, because due to the RTL hack to cut
           text off on left side the punctuation characters get
           reordered (of course). Prepending &lrm; works for most
           characters, parentheses we need to *soak* in it. But only
           the right ones. And that for some reason needs to be also for &.
           Huh. https://en.wikipedia.org/wiki/Right-to-left_mark */
        return this.escape(name).replace(/[:=]/g, '&lrm;$&').replace(/(\)|&gt;|&amp;|\/)/g, '&lrm;$&&lrm;');
    },

    renderResults: /* istanbul ignore next */ function(resultsSuggestedTabAutocompletion) {
        if(!this.searchString.length) {
            document.getElementById('search-help').style.display = 'block';
            document.getElementById('search-results').style.display = 'none';
            document.getElementById('search-notfound').style.display = 'none';
            return;
        }

        document.getElementById('search-help').style.display = 'none';

        /* Results found */
        if(resultsSuggestedTabAutocompletion[0].length) {
            let results = resultsSuggestedTabAutocompletion[0];

            document.getElementById('search-results').style.display = 'block';
            document.getElementById('search-notfound').style.display = 'none';

            let list = '';
            for(let i = 0; i != results.length; ++i) {
                /* Labels + */
                list += '<li' + (i ? '' : ' id="search-current"') + '><a href="' + results[i].url + '" onmouseover="selectResult(event)" data-md-link-title="' + this.escape(results[i].name.substr(results[i].name.length - this.searchString.length - results[i].suffixLength)) + '"><div class="m-label m-flat ' + results[i].cssClass + '">' + results[i].typeName + '</div>' + (results[i].flags & 2 ? '<div class="m-label m-danger">deprecated</div>' : '') + (results[i].flags & 4 ? '<div class="m-label m-danger">deleted</div>' : '');

                /* Render the alias (cut off from the right) */
                if(results[i].alias) {
                    list += '<div class="m-doc-search-alias"><span class="m-text m-dim">' + this.escape(results[i].name.substr(0, results[i].name.length - this.searchString.length - results[i].suffixLength)) + '</span><span class="m-doc-search-typed">' + this.escape(results[i].name.substr(results[i].name.length - this.searchString.length - results[i].suffixLength, this.searchString.length)) + '</span>' + this.escapeForRtl(results[i].name.substr(results[i].name.length - results[i].suffixLength)) + '<span class="m-text m-dim">: ' + this.escape(results[i].alias) + '</span>';

                /* Render the normal thing (cut off from the left, have to
                   escape for RTL) */
                } else {
                    list += '<div><span class="m-text m-dim">' + this.escapeForRtl(results[i].name.substr(0, results[i].name.length - this.searchString.length - results[i].suffixLength)) + '</span><span class="m-doc-search-typed">' + this.escapeForRtl(results[i].name.substr(results[i].name.length - this.searchString.length - results[i].suffixLength, this.searchString.length)) + '</span>' + this.escapeForRtl(results[i].name.substr(results[i].name.length - results[i].suffixLength));
                }

                /* The closing */
                list += '</div></a></li>';
            }
            document.getElementById('search-results').innerHTML = this.fromUtf8(list);
            document.getElementById('search-current').scrollIntoView(true);

            /* Append the suggested tab autocompletion, if any, and if the user
               didn't just delete it */
            let searchInput = document.getElementById('search-input');
            if(this.autocompleteNextInputEvent && resultsSuggestedTabAutocompletion[1].length && searchInput.selectionEnd == searchInput.value.length) {
                let suggestedTabAutocompletion = this.fromUtf8(resultsSuggestedTabAutocompletion[1]);

                let lengthBefore = searchInput.value.length;
                searchInput.value += suggestedTabAutocompletion;
                searchInput.setSelectionRange(lengthBefore, searchInput.value.length);
            }

        /* Nothing found */
        } else {
            document.getElementById('search-results').innerHTML = '';
            document.getElementById('search-results').style.display = 'none';
            document.getElementById('search-notfound').style.display = 'block';
        }

        /* Don't allow things to be selected just by motionless mouse cursor
           suddenly appearing over a search result */
        this.mouseMovedSinceLastRender = false;

        /* Reset autocompletion, if it was allowed. It'll get whitelisted next
           time a character gets inserted. */
        this.autocompleteNextInputEvent = false;
    },

    searchAndRender: /* istanbul ignore next */ function(value) {
        let prev = performance.now();
        let results = this.search(value);
        let after = performance.now();
        this.renderResults(results);
        if(this.searchString.length) {
            document.getElementById('search-symbolcount').innerHTML =
                results[0].length + (results[0].length >= this.maxResults ? '+' : '') + " results (" + Math.round((after - prev)*10)/10 + " ms)";
        } else
            document.getElementById('search-symbolcount').innerHTML = this.symbolCount;
    },
};

/* istanbul ignore next */
function selectResult(event) {
    if(!Search.mouseMovedSinceLastRender) return;

    if(event.currentTarget.parentNode.id == 'search-current') return;

    let current = document.getElementById('search-current');
    current.removeAttribute('id');
    event.currentTarget.parentNode.id = 'search-current';
}

/* This is separated from showSearch() because we need non-destructive behavior
   when appearing directly on a URL with #search */ /* istanbul ignore next */
function updateForSearchVisible() {
    /* Prevent accidental scrolling of the body, prevent page layout jumps */
    let scrolledBodyWidth = document.body.offsetWidth;
    document.body.style.overflow = 'hidden';
    document.body.style.paddingRight = (document.body.offsetWidth - scrolledBodyWidth) + 'px';

    document.getElementById('search-input').value = '';
    document.getElementById('search-input').focus();
    document.getElementById('search-results').style.display = 'none';
    document.getElementById('search-notfound').style.display = 'none';
    document.getElementById('search-help').style.display = 'block';
}

/* istanbul ignore next */
function showSearch() {
    window.location.hash = '#search';
    Search.canGoBackToHideSearch = true;

    updateForSearchVisible();
    document.getElementById('search-symbolcount').innerHTML = Search.symbolCount;
    return false;
}

/* istanbul ignore next */
function hideSearch() {
    /* If the search box was opened using showSearch(), we can go back in the
       history. Otherwise (for example when we landed to #search from a
       bookmark or another server), going back would not do the right thing and
       in that case we simply replace the current history state. */
    if(Search.canGoBackToHideSearch) {
        Search.canGoBackToHideSearch = false;
        window.history.back();
    } else {
        window.location.hash = '#!';
        window.history.replaceState('', '', window.location.pathname);
    }

    /* Restore scrollbar, prevent page layout jumps */
    document.body.style.overflow = 'auto';
    document.body.style.paddingRight = '0';

    return false;
}

/* istanbul ignore next */
function copyToKeyboard(text) {
    /* Append to the popup, appending to document.body would cause it to
       scroll when focused */
    let searchPopup = document.getElementsByClassName('m-doc-search')[0];
    let textarea = document.createElement("textarea");
    textarea.value = text;
    searchPopup.appendChild(textarea);
    textarea.focus();
    textarea.select();

    document.execCommand('copy');

    searchPopup.removeChild(textarea);
    document.getElementById('search-input').focus();
}

/* Only in case we're running in a browser. Why a simple if(document) doesn't
   work is beyond me. */ /* istanbul ignore if */
if(typeof document !== 'undefined') {
    document.getElementById('search-input').oninput = function(event) {
        Search.searchAndRender(document.getElementById('search-input').value);
    };

    document.onkeydown = function(event) {
        /* Search shown */
        if(window.location.hash == '#search') {
            /* Close the search */
            if(event.key == 'Escape') {
                hideSearch();

            /* Focus the search input, if not already, using T or Tab */
            } else if((!document.activeElement || document.activeElement.id != 'search-input') && (event.key.toLowerCase() == 't' || event.key == 'Tab') && !event.shiftKey && !event.ctrlKey && !event.altKey && !event.metaKey) {
                document.getElementById('search-input').focus();
                return false; /* so T doesn't get entered into the box */

            /* Fill in the autocompleted selection */
            } else if(event.key == 'Tab' && !event.shiftKey && !event.ctrlKey && !event.altKey && !event.metaKey) {
                /* But only if the input has selection at the end */
                let input = document.getElementById('search-input');
                if(input.selectionEnd == input.value.length && input.selectionStart != input.selectionEnd) {
                    input.setSelectionRange(input.value.length, input.value.length);
                    return false; /* so input won't lose focus */
                }

            /* Select next item */
            } else if(event.key == 'ArrowDown') {
                let current = document.getElementById('search-current');
                if(current) {
                    let next = current.nextSibling;
                    if(next) {
                        current.id = '';
                        next.id = 'search-current';
                        next.scrollIntoView(false);
                    }
                }
                return false; /* so the keypress doesn't affect input cursor */

            /* Select prev item */
            } else if(event.key == 'ArrowUp') {
                let current = document.getElementById('search-current');
                if(current) {
                    let prev = current.previousSibling;
                    if(prev) {
                        current.id = '';
                        prev.id = 'search-current';
                        prev.scrollIntoView(false);
                    }
                }
                return false; /* so the keypress doesn't affect input cursor */

            /* Go to result (if any) */
            } else if(event.key == 'Enter') {
                let result = document.getElementById('search-current');
                if(result) {
                    result.firstElementChild.click();

                    /* We might be staying on the same page, so restore scrollbar,
                       and prevent page layout jumps */
                    document.body.style.overflow = 'auto';
                    document.body.style.paddingRight = '0';
                }
                return false; /* so the form doesn't get sent */

            /* Copy (Markdown) link to keyboard */
            } else if((event.key.toLowerCase() == 'l' || event.key.toLowerCase() == 'm') && event.metaKey) {
                let result = document.getElementById('search-current');
                if(result) {
                    let plain = event.key.toLowerCase() == 'l';
                    let link = plain ? result.firstElementChild.href :
                        '[' + result.firstElementChild.dataset.mdLinkTitle + '](' + result.firstElementChild.href + ')';

                    copyToKeyboard(link);

                    /* Add CSS class to the element for visual feedback (this
                       will get removed on keyup), but only if it's not already
                       there (in case of key repeat, e.g.) */
                    if(result.className.indexOf('m-doc-search-copied') == -1)
                        result.className += ' m-doc-search-copied';
                    console.log("Copied " +  (plain ? "link" : "Markdown link") + " to " + result.firstElementChild.dataset.mdLinkTitle);
                }

                return false; /* so L doesn't get entered into the box */

            /* Looks like the user is inserting some text (and not cutting,
               copying or whatever), allow autocompletion for the new
               character. The oninput event resets this back to false, so this
               basically whitelists only keyboard input, including Shift-key
               and special chars using right Alt (or equivalent on Mac), but
               excluding Ctrl-key, which is usually not for text input. In the
               worst case the autocompletion won't be allowed ever, which is
               much more acceptable behavior than having no ability to disable
               it and annoying the users. */
            } else if(event.key != 'Backspace' && event.key != 'Delete' && !event.metaKey && (!event.ctrlKey || event.altKey)
                /* Don't ever attempt autocompletion with Android virtual
                   keyboards, as those report all `event.key`s as
                   `Unidentified` (on Chrome) or `Process` (on Firefox) with
                   `event.code` 229 and thus we have no way to tell if a text
                   is entered or deleted. See this WONTFIX bug for details:
                    https://bugs.chromium.org/p/chromium/issues/detail?id=118639
                   Couldn't find any similar bugreport for Firefox, but I
                   assume the virtual keyboard is to blame.

                   An alternative is to hook into inputEvent, which has the
                   data, but ... there's more cursed issues right after that:

                    - setSelectionRange() in Chrome on Android only renders
                      stuff, but doesn't actually act as such. Pressing
                      Backspace will only remove the highlight, but the text
                      stays here. Only delay-calling it through a timeout will
                      work as intended. Possibly related SO suggestion (back
                      then not even the rendering worked properly):
                       https://stackoverflow.com/a/13235951
                      Possibly related Chrome bug:
                       https://bugs.chromium.org/p/chromium/issues/detail?id=32865

                    - On Firefox Mobile, programmatically changing an input
                      value (for the autocompletion highlight) will trigger an
                      input event, leading to search *and* autocompletion being
                      triggered again. Ultimately that results in newly typed
                      characters not replacing the autocompletion but rather
                      inserting before it, corrupting the searched string. This
                      event has to be explicitly ignored.

                    - On Firefox Mobile, deleting a highlight with the
                      backspace key will result in *three* input events instead
                      of one:
                        1. `deleteContentBackward` removing the selection (same
                           as Chrome or desktop Firefox)
                        2. `deleteContentBackward` removing *the whole word*
                           that contained the selection (or the whole text if
                           it's just one word)
                        3. `insertCompositionText`, adding the word back in,
                           resulting in the same state as (1).
                      I have no idea WHY it has to do this (possibly some
                      REALLY NASTY workaround to trigger correct font shaping?)
                      but ultimately it results in the autocompletion being
                      added again right after it got deleted, making this whole
                      thing VERY annoying to use.

                   I attempted to work around the above, but it resulted in a
                   huge amount of browser-specific code that achieves only 90%
                   of the goal, with certain corner cases still being rather
                   broken (such as autocompletion randomly triggering when
                   erasing the text, even though it shouldn't). So disabling
                   autocompletion on this HELLISH BROKEN PLATFORM is the best
                   option at the moment. */
                && event.key != 'Unidentified' && event.key != 'Process'
            ) {
                Search.autocompleteNextInputEvent = true;
            /* Otherwise reset the flag, because when the user would press e.g.
               the 'a' key and then e.g. ArrowRight (which doesn't trigger
               oninput), a Backspace after would still result in
               autocompleteNextInputEvent, because nothing reset it back. */
            } else {
                Search.autocompleteNextInputEvent = false;
            }

        /* Search hidden */
        } else {
            /* Open the search on the T or Tab key */
            if((event.key.toLowerCase() == 't' || event.key == 'Tab') && !event.shiftKey && !event.ctrlKey && !event.altKey && !event.metaKey) {
                showSearch();
                return false; /* so T doesn't get entered into the box */
            }
        }
    };

    document.onkeyup = function(event) {
        /* Remove highlight after key is released after a link copy */
        if((event.key.toLowerCase() == 'l' || event.key.toLowerCase() == 'm') && event.metaKey) {
            let result = document.getElementById('search-current');
            if(result) result.className = result.className.replace(' m-doc-search-copied', '');
        }
    };

    /* Allow selecting items by mouse hover only after it moves once the
       results are populated. This prevents a random item getting selected if
       the cursor is left motionless over the result area. */
    document.getElementById('search-results').onmousemove = function() {
        Search.mouseMovedSinceLastRender = true;
    };

    /* If #search is already present in the URL, hide the scrollbar etc. for a
       consistent experience */
    if(window.location.hash == '#search') updateForSearchVisible();
}

/* For Node.js testing */ /* istanbul ignore else */
if(typeof module !== 'undefined') { module.exports = { Search: Search }; }

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
})(Search, [{"name": "Rodin::Adaptation::CellDeformation", "url": "class_rodin_1_1_adaptation_1_1_cell_deformation.html", "keys": ["rodin::adaptation::celldeformation", "adaptation::celldeformation", "celldeformation"]}, {"name": "Rodin::Adaptation::DeformationMap", "url": "class_rodin_1_1_adaptation_1_1_deformation_map.html", "keys": ["rodin::adaptation::deformationmap", "adaptation::deformationmap", "deformationmap"]}, {"name": "Rodin::Adaptation::SWIFT::Adapt", "url": "class_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_adapt.html", "keys": ["rodin::adaptation::swift::adapt", "adaptation::swift::adapt", "swift::adapt", "adapt"]}, {"name": "Rodin::Adaptation::SWIFT::CurrentInverse", "url": "class_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_current_inverse.html", "keys": ["rodin::adaptation::swift::currentinverse", "adaptation::swift::currentinverse", "swift::currentinverse", "currentinverse"]}, {"name": "Rodin::Adaptation::SWIFT::Distribution", "url": "class_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_distribution.html", "keys": ["rodin::adaptation::swift::distribution", "adaptation::swift::distribution", "swift::distribution", "distribution"]}, {"name": "Rodin::Adaptation::SWIFT::FittingForce", "url": "class_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_fitting_force.html", "keys": ["rodin::adaptation::swift::fittingforce", "adaptation::swift::fittingforce", "swift::fittingforce", "fittingforce"]}, {"name": "Rodin::Adaptation::SWIFT::FittingTensor", "url": "class_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_fitting_tensor.html", "keys": ["rodin::adaptation::swift::fittingtensor", "adaptation::swift::fittingtensor", "swift::fittingtensor", "fittingtensor"]}, {"name": "Rodin::Adaptation::SWIFT::HingeForce", "url": "class_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_hinge_force.html", "keys": ["rodin::adaptation::swift::hingeforce", "adaptation::swift::hingeforce", "swift::hingeforce", "hingeforce"]}, {"name": "Rodin::Adaptation::SWIFT::HingeMetric", "url": "class_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_hinge_metric.html", "keys": ["rodin::adaptation::swift::hingemetric", "adaptation::swift::hingemetric", "swift::hingemetric", "hingemetric"]}, {"name": "Rodin::Adaptation::SWIFT::HingeProblem", "url": "class_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_hinge_problem.html", "keys": ["rodin::adaptation::swift::hingeproblem", "adaptation::swift::hingeproblem", "swift::hingeproblem", "hingeproblem"]}, {"name": "Rodin::Adaptation::SWIFT::HingeState", "url": "class_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_hinge_state.html", "keys": ["rodin::adaptation::swift::hingestate", "adaptation::swift::hingestate", "swift::hingestate", "hingestate"]}, {"name": "Rodin::Adaptation::SWIFT::LinearSolver", "url": "class_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_linear_solver.html", "keys": ["rodin::adaptation::swift::linearsolver", "adaptation::swift::linearsolver", "swift::linearsolver", "linearsolver"]}, {"name": "Rodin::Adaptation::SWIFT::Loss", "url": "class_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_loss.html", "keys": ["rodin::adaptation::swift::loss", "adaptation::swift::loss", "swift::loss", "loss"]}, {"name": "Rodin::Adaptation::SWIFT::Problem", "url": "class_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_problem.html", "keys": ["rodin::adaptation::swift::problem", "adaptation::swift::problem", "swift::problem", "problem"]}, {"name": "Rodin::Adaptation::SWIFT::QualityCovering", "url": "class_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_quality_covering.html", "keys": ["rodin::adaptation::swift::qualitycovering", "adaptation::swift::qualitycovering", "swift::qualitycovering", "qualitycovering"]}, {"name": "Rodin::Adaptation::SWIFT::QualitySamples", "url": "class_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_quality_samples.html", "keys": ["rodin::adaptation::swift::qualitysamples", "adaptation::swift::qualitysamples", "swift::qualitysamples", "qualitysamples"]}, {"name": "Rodin::Adaptation::SWIFT::ResidualState", "url": "class_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_residual_state.html", "keys": ["rodin::adaptation::swift::residualstate", "adaptation::swift::residualstate", "swift::residualstate", "residualstate"]}, {"name": "Rodin::Adaptation::SWIFT::WitnessLattice", "url": "class_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_witness_lattice.html", "keys": ["rodin::adaptation::swift::witnesslattice", "adaptation::swift::witnesslattice", "swift::witnesslattice", "witnesslattice"]}, {"name": "Rodin::Advection::Lagrangian", "url": "class_rodin_1_1_advection_1_1_lagrangian.html", "keys": ["rodin::advection::lagrangian", "advection::lagrangian", "lagrangian"]}, {"name": "Rodin::Advection::StopInsideBoundaryPolicy", "url": "class_rodin_1_1_advection_1_1_stop_inside_boundary_policy.html", "keys": ["rodin::advection::stopinsideboundarypolicy", "advection::stopinsideboundarypolicy", "stopinsideboundarypolicy"]}, {"name": "Rodin::Advection::TaylorBoundaryShiftPolicy", "url": "class_rodin_1_1_advection_1_1_taylor_boundary_shift_policy.html", "keys": ["rodin::advection::taylorboundaryshiftpolicy", "advection::taylorboundaryshiftpolicy", "taylorboundaryshiftpolicy"]}, {"name": "Rodin::Alert::ClassException", "url": "class_rodin_1_1_alert_1_1_class_exception.html", "keys": ["rodin::alert::classexception", "alert::classexception", "classexception"]}, {"name": "Rodin::Alert::Color", "url": "class_rodin_1_1_alert_1_1_color.html", "keys": ["rodin::alert::color", "alert::color", "color"]}, {"name": "Rodin::Alert::Exception", "url": "class_rodin_1_1_alert_1_1_exception.html", "keys": ["rodin::alert::exception", "alert::exception", "exception"]}, {"name": "Rodin::Alert::ExceptionPrefix", "url": "class_rodin_1_1_alert_1_1_exception_prefix.html", "keys": ["rodin::alert::exceptionprefix", "alert::exceptionprefix", "exceptionprefix"]}, {"name": "Rodin::Alert::IdentifierT", "url": "class_rodin_1_1_alert_1_1_identifier_t.html", "keys": ["rodin::alert::identifiert", "alert::identifiert", "identifiert"]}, {"name": "Rodin::Alert::Info", "url": "class_rodin_1_1_alert_1_1_info.html", "keys": ["rodin::alert::info", "alert::info", "info"]}, {"name": "Rodin::Alert::InfoPrefix", "url": "class_rodin_1_1_alert_1_1_info_prefix.html", "keys": ["rodin::alert::infoprefix", "alert::infoprefix", "infoprefix"]}, {"name": "Rodin::Alert::Internal::CanBeOutput", "url": "class_rodin_1_1_alert_1_1_internal_1_1_can_be_output.html", "keys": ["rodin::alert::internal::canbeoutput", "alert::internal::canbeoutput", "internal::canbeoutput", "canbeoutput"]}, {"name": "Rodin::Alert::MemberFunctionException", "url": "class_rodin_1_1_alert_1_1_member_function_exception.html", "keys": ["rodin::alert::memberfunctionexception", "alert::memberfunctionexception", "memberfunctionexception"]}, {"name": "Rodin::Alert::MemberFunctionWarning", "url": "class_rodin_1_1_alert_1_1_member_function_warning.html", "keys": ["rodin::alert::memberfunctionwarning", "alert::memberfunctionwarning", "memberfunctionwarning"]}, {"name": "Rodin::Alert::Message", "url": "class_rodin_1_1_alert_1_1_message.html", "keys": ["rodin::alert::message", "alert::message", "message"]}, {"name": "Rodin::Alert::MessagePrefix", "url": "class_rodin_1_1_alert_1_1_message_prefix.html", "keys": ["rodin::alert::messageprefix", "alert::messageprefix", "messageprefix"]}, {"name": "Rodin::Alert::NamespacedException", "url": "class_rodin_1_1_alert_1_1_namespaced_exception.html", "keys": ["rodin::alert::namespacedexception", "alert::namespacedexception", "namespacedexception"]}, {"name": "Rodin::Alert::Notation", "url": "class_rodin_1_1_alert_1_1_notation.html", "keys": ["rodin::alert::notation", "alert::notation", "notation"]}, {"name": "Rodin::Alert::PrefixedMessage", "url": "class_rodin_1_1_alert_1_1_prefixed_message.html", "keys": ["rodin::alert::prefixedmessage", "alert::prefixedmessage", "prefixedmessage"]}, {"name": "Rodin::Alert::Success", "url": "class_rodin_1_1_alert_1_1_success.html", "keys": ["rodin::alert::success", "alert::success", "success"]}, {"name": "Rodin::Alert::SuccessPrefix", "url": "class_rodin_1_1_alert_1_1_success_prefix.html", "keys": ["rodin::alert::successprefix", "alert::successprefix", "successprefix"]}, {"name": "Rodin::Alert::Text", "url": "class_rodin_1_1_alert_1_1_text.html", "keys": ["rodin::alert::text", "alert::text", "text"]}, {"name": "Rodin::Alert::Warning", "url": "class_rodin_1_1_alert_1_1_warning.html", "keys": ["rodin::alert::warning", "alert::warning", "warning"]}, {"name": "Rodin::Alert::WarningPrefix", "url": "class_rodin_1_1_alert_1_1_warning_prefix.html", "keys": ["rodin::alert::warningprefix", "alert::warningprefix", "warningprefix"]}, {"name": "Rodin::Assembly::AssemblyBase", "url": "class_rodin_1_1_assembly_1_1_assembly_base.html", "keys": ["rodin::assembly::assemblybase", "assembly::assemblybase", "assemblybase"]}, {"name": "Rodin::Assembly::BilinearFormAssemblyInput", "url": "class_rodin_1_1_assembly_1_1_bilinear_form_assembly_input.html", "keys": ["rodin::assembly::bilinearformassemblyinput", "assembly::bilinearformassemblyinput", "bilinearformassemblyinput"]}, {"name": "Rodin::Assembly::BilinearFormTupleAssemblyInput", "url": "class_rodin_1_1_assembly_1_1_bilinear_form_tuple_assembly_input.html", "keys": ["rodin::assembly::bilinearformtupleassemblyinput", "assembly::bilinearformtupleassemblyinput", "bilinearformtupleassemblyinput"]}, {"name": "Rodin::Assembly::ConstraintMap", "url": "class_rodin_1_1_assembly_1_1_constraint_map.html", "keys": ["rodin::assembly::constraintmap", "assembly::constraintmap", "constraintmap"]}, {"name": "Rodin::Assembly::Default", "url": "class_rodin_1_1_assembly_1_1_default.html", "keys": ["rodin::assembly::default", "assembly::default", "default"]}, {"name": "Rodin::Assembly::DirichletBCAssemblyInput", "url": "class_rodin_1_1_assembly_1_1_dirichlet_b_c_assembly_input.html", "keys": ["rodin::assembly::dirichletbcassemblyinput", "assembly::dirichletbcassemblyinput", "dirichletbcassemblyinput"]}, {"name": "Rodin::Assembly::DirichletBCShapeFunctionAssemblyInput", "url": "class_rodin_1_1_assembly_1_1_dirichlet_b_c_shape_function_assembly_input.html", "keys": ["rodin::assembly::dirichletbcshapefunctionassemblyinput", "assembly::dirichletbcshapefunctionassemblyinput", "dirichletbcshapefunctionassemblyinput"]}, {"name": "Rodin::Assembly::LinearFormAssemblyInput", "url": "class_rodin_1_1_assembly_1_1_linear_form_assembly_input.html", "keys": ["rodin::assembly::linearformassemblyinput", "assembly::linearformassemblyinput", "linearformassemblyinput"]}, {"name": "Rodin::Assembly::LinearFormTupleAssemblyInput", "url": "class_rodin_1_1_assembly_1_1_linear_form_tuple_assembly_input.html", "keys": ["rodin::assembly::linearformtupleassemblyinput", "assembly::linearformtupleassemblyinput", "linearformtupleassemblyinput"]}, {"name": "Rodin::Assembly::MPI", "url": "class_rodin_1_1_assembly_1_1_m_p_i.html", "keys": ["rodin::assembly::mpi", "assembly::mpi", "mpi"]}, {"name": "Rodin::Assembly::MPIBoundaryDOFs", "url": "class_rodin_1_1_assembly_1_1_m_p_i_boundary_d_o_fs.html", "keys": ["rodin::assembly::mpiboundarydofs", "assembly::mpiboundarydofs", "mpiboundarydofs"]}, {"name": "Rodin::Assembly::MPIIteration", "url": "class_rodin_1_1_assembly_1_1_m_p_i_iteration.html", "keys": ["rodin::assembly::mpiiteration", "assembly::mpiiteration", "mpiiteration"]}, {"name": "Rodin::Assembly::OpenMP", "url": "class_rodin_1_1_assembly_1_1_open_m_p.html", "keys": ["rodin::assembly::openmp", "assembly::openmp", "openmp"]}, {"name": "Rodin::Assembly::OpenMPIteration", "url": "class_rodin_1_1_assembly_1_1_open_m_p_iteration.html", "keys": ["rodin::assembly::openmpiteration", "assembly::openmpiteration", "openmpiteration"]}, {"name": "Rodin::Assembly::ProblemAssemblyInput", "url": "class_rodin_1_1_assembly_1_1_problem_assembly_input.html", "keys": ["rodin::assembly::problemassemblyinput", "assembly::problemassemblyinput", "problemassemblyinput"]}, {"name": "Rodin::Assembly::ScatterMap", "url": "class_rodin_1_1_assembly_1_1_scatter_map.html", "keys": ["rodin::assembly::scattermap", "assembly::scattermap", "scattermap"]}, {"name": "Rodin::Assembly::Sequential", "url": "class_rodin_1_1_assembly_1_1_sequential.html", "keys": ["rodin::assembly::sequential", "assembly::sequential", "sequential"]}, {"name": "Rodin::Assembly::SequentialIteration", "url": "class_rodin_1_1_assembly_1_1_sequential_iteration.html", "keys": ["rodin::assembly::sequentialiteration", "assembly::sequentialiteration", "sequentialiteration"]}, {"name": "Rodin::Cast", "url": "class_rodin_1_1_cast.html", "keys": ["rodin::cast", "cast"]}, {"name": "Rodin::Context::Base", "url": "class_rodin_1_1_context_1_1_base.html", "keys": ["rodin::context::base", "context::base", "base"]}, {"name": "Rodin::Context::Local", "url": "class_rodin_1_1_context_1_1_local.html", "keys": ["rodin::context::local", "context::local", "local"]}, {"name": "Rodin::Context::MPI", "url": "class_rodin_1_1_context_1_1_m_p_i.html", "keys": ["rodin::context::mpi", "context::mpi", "mpi"]}, {"name": "Rodin::Copyable", "url": "class_rodin_1_1_copyable.html", "keys": ["rodin::copyable", "copyable"]}, {"name": "Rodin::Distance::Base", "url": "class_rodin_1_1_distance_1_1_base.html", "keys": ["rodin::distance::base", "distance::base", "base"]}, {"name": "Rodin::Distance::Eikonal", "url": "class_rodin_1_1_distance_1_1_eikonal.html", "keys": ["rodin::distance::eikonal", "distance::eikonal", "eikonal"]}, {"name": "Rodin::Distance::Poisson", "url": "class_rodin_1_1_distance_1_1_poisson.html", "keys": ["rodin::distance::poisson", "distance::poisson", "poisson"]}, {"name": "Rodin::Distance::Rvachev", "url": "class_rodin_1_1_distance_1_1_rvachev.html", "keys": ["rodin::distance::rvachev", "distance::rvachev", "rvachev"]}, {"name": "Rodin::Distance::SignedPoisson", "url": "class_rodin_1_1_distance_1_1_signed_poisson.html", "keys": ["rodin::distance::signedpoisson", "distance::signedpoisson", "signedpoisson"]}, {"name": "Rodin::Distance::SpaldingTucker", "url": "class_rodin_1_1_distance_1_1_spalding_tucker.html", "keys": ["rodin::distance::spaldingtucker", "distance::spaldingtucker", "spaldingtucker"]}, {"name": "Rodin::Eikonal::FMM", "url": "class_rodin_1_1_eikonal_1_1_f_m_m.html", "keys": ["rodin::eikonal::fmm", "eikonal::fmm", "fmm"]}, {"name": "Rodin::FormLanguage::Base", "url": "class_rodin_1_1_form_language_1_1_base.html", "keys": ["rodin::formlanguage::base", "formlanguage::base", "base"]}, {"name": "Rodin::FormLanguage::List", "url": "class_rodin_1_1_form_language_1_1_list.html", "keys": ["rodin::formlanguage::list", "formlanguage::list", "list"]}, {"name": "Rodin::FormLanguage::List::ConstIterator", "url": "class_rodin_1_1_form_language_1_1_list_1_1_const_iterator.html", "keys": ["rodin::formlanguage::list::constiterator", "formlanguage::list::constiterator", "list::constiterator", "constiterator"]}, {"name": "Rodin::FormLanguage::List::Iterator", "url": "class_rodin_1_1_form_language_1_1_list_1_1_iterator.html", "keys": ["rodin::formlanguage::list::iterator", "formlanguage::list::iterator", "list::iterator", "iterator"]}, {"name": "Rodin::Geometry::AttributeIndex", "url": "class_rodin_1_1_geometry_1_1_attribute_index.html", "keys": ["rodin::geometry::attributeindex", "geometry::attributeindex", "attributeindex"]}, {"name": "Rodin::Geometry::BalancedCompactPartitioner", "url": "class_rodin_1_1_geometry_1_1_balanced_compact_partitioner.html", "keys": ["rodin::geometry::balancedcompactpartitioner", "geometry::balancedcompactpartitioner", "balancedcompactpartitioner"]}, {"name": "Rodin::Geometry::BoundedIndexGenerator", "url": "class_rodin_1_1_geometry_1_1_bounded_index_generator.html", "keys": ["rodin::geometry::boundedindexgenerator", "geometry::boundedindexgenerator", "boundedindexgenerator"]}, {"name": "Rodin::Geometry::Bounds", "url": "class_rodin_1_1_geometry_1_1_bounds.html", "keys": ["rodin::geometry::bounds", "geometry::bounds", "bounds"]}, {"name": "Rodin::Geometry::CCL", "url": "class_rodin_1_1_geometry_1_1_c_c_l.html", "keys": ["rodin::geometry::ccl", "geometry::ccl", "ccl"]}, {"name": "Rodin::Geometry::Cell", "url": "class_rodin_1_1_geometry_1_1_cell.html", "keys": ["rodin::geometry::cell", "geometry::cell", "cell"]}, {"name": "Rodin::Geometry::CellIterator", "url": "class_rodin_1_1_geometry_1_1_cell_iterator.html", "keys": ["rodin::geometry::celliterator", "geometry::celliterator", "celliterator"]}, {"name": "Rodin::Geometry::Connectivity", "url": "class_rodin_1_1_geometry_1_1_connectivity.html", "keys": ["rodin::geometry::connectivity", "geometry::connectivity", "connectivity"]}, {"name": "Rodin::Geometry::ConnectivityBase", "url": "class_rodin_1_1_geometry_1_1_connectivity_base.html", "keys": ["rodin::geometry::connectivitybase", "geometry::connectivitybase", "connectivitybase"]}, {"name": "Rodin::Geometry::ContainerIndexGenerator", "url": "class_rodin_1_1_geometry_1_1_container_index_generator.html", "keys": ["rodin::geometry::containerindexgenerator", "geometry::containerindexgenerator", "containerindexgenerator"]}, {"name": "Rodin::Geometry::EmptyIndexGenerator", "url": "class_rodin_1_1_geometry_1_1_empty_index_generator.html", "keys": ["rodin::geometry::emptyindexgenerator", "geometry::emptyindexgenerator", "emptyindexgenerator"]}, {"name": "Rodin::Geometry::Euclidean::Base", "url": "class_rodin_1_1_geometry_1_1_euclidean_1_1_base.html", "keys": ["rodin::geometry::euclidean::base", "geometry::euclidean::base", "euclidean::base", "base"]}, {"name": "Rodin::Geometry::Euclidean::Box", "url": "class_rodin_1_1_geometry_1_1_euclidean_1_1_box.html", "keys": ["rodin::geometry::euclidean::box", "geometry::euclidean::box", "euclidean::box", "box"]}, {"name": "Rodin::Geometry::Euclidean::Hyperplane", "url": "class_rodin_1_1_geometry_1_1_euclidean_1_1_hyperplane.html", "keys": ["rodin::geometry::euclidean::hyperplane", "geometry::euclidean::hyperplane", "euclidean::hyperplane", "hyperplane"]}, {"name": "Rodin::Geometry::Euclidean::Line", "url": "class_rodin_1_1_geometry_1_1_euclidean_1_1_line.html", "keys": ["rodin::geometry::euclidean::line", "geometry::euclidean::line", "euclidean::line", "line"]}, {"name": "Rodin::Geometry::Euclidean::LineSegment", "url": "class_rodin_1_1_geometry_1_1_euclidean_1_1_line_segment.html", "keys": ["rodin::geometry::euclidean::linesegment", "geometry::euclidean::linesegment", "euclidean::linesegment", "linesegment"]}, {"name": "Rodin::Geometry::Euclidean::Point", "url": "class_rodin_1_1_geometry_1_1_euclidean_1_1_point.html", "keys": ["rodin::geometry::euclidean::point", "geometry::euclidean::point", "euclidean::point", "point"]}, {"name": "Rodin::Geometry::Euclidean::Sphere", "url": "class_rodin_1_1_geometry_1_1_euclidean_1_1_sphere.html", "keys": ["rodin::geometry::euclidean::sphere", "geometry::euclidean::sphere", "euclidean::sphere", "sphere"]}, {"name": "Rodin::Geometry::Face", "url": "class_rodin_1_1_geometry_1_1_face.html", "keys": ["rodin::geometry::face", "geometry::face", "face"]}, {"name": "Rodin::Geometry::FaceIterator", "url": "class_rodin_1_1_geometry_1_1_face_iterator.html", "keys": ["rodin::geometry::faceiterator", "geometry::faceiterator", "faceiterator"]}, {"name": "Rodin::Geometry::GeometryIndexed", "url": "class_rodin_1_1_geometry_1_1_geometry_indexed.html", "keys": ["rodin::geometry::geometryindexed", "geometry::geometryindexed", "geometryindexed"]}, {"name": "Rodin::Geometry::GreedyPartitioner", "url": "class_rodin_1_1_geometry_1_1_greedy_partitioner.html", "keys": ["rodin::geometry::greedypartitioner", "geometry::greedypartitioner", "greedypartitioner"]}, {"name": "Rodin::Geometry::IdentityTransformation", "url": "class_rodin_1_1_geometry_1_1_identity_transformation.html", "keys": ["rodin::geometry::identitytransformation", "geometry::identitytransformation", "identitytransformation"]}, {"name": "Rodin::Geometry::IndexGeneratorBase", "url": "class_rodin_1_1_geometry_1_1_index_generator_base.html", "keys": ["rodin::geometry::indexgeneratorbase", "geometry::indexgeneratorbase", "indexgeneratorbase"]}, {"name": "Rodin::Geometry::IteratorIndexGenerator", "url": "class_rodin_1_1_geometry_1_1_iterator_index_generator.html", "keys": ["rodin::geometry::iteratorindexgenerator", "geometry::iteratorindexgenerator", "iteratorindexgenerator"]}, {"name": "Rodin::Geometry::Mesh", "url": "class_rodin_1_1_geometry_1_1_mesh.html", "keys": ["rodin::geometry::mesh", "geometry::mesh", "mesh"]}, {"name": "Rodin::Geometry::MeshBase", "url": "class_rodin_1_1_geometry_1_1_mesh_base.html", "keys": ["rodin::geometry::meshbase", "geometry::meshbase", "meshbase"]}, {"name": "Rodin::Geometry::MinSTCut", "url": "class_rodin_1_1_geometry_1_1_min_s_t_cut.html", "keys": ["rodin::geometry::minstcut", "geometry::minstcut", "minstcut"]}, {"name": "Rodin::Geometry::ParametricTransformation", "url": "class_rodin_1_1_geometry_1_1_parametric_transformation.html", "keys": ["rodin::geometry::parametrictransformation", "geometry::parametrictransformation", "parametrictransformation"]}, {"name": "Rodin::Geometry::Partitioner", "url": "class_rodin_1_1_geometry_1_1_partitioner.html", "keys": ["rodin::geometry::partitioner", "geometry::partitioner", "partitioner"]}, {"name": "Rodin::Geometry::Point", "url": "class_rodin_1_1_geometry_1_1_point.html", "keys": ["rodin::geometry::point", "geometry::point", "point"]}, {"name": "Rodin::Geometry::PointBase", "url": "class_rodin_1_1_geometry_1_1_point_base.html", "keys": ["rodin::geometry::pointbase", "geometry::pointbase", "pointbase"]}, {"name": "Rodin::Geometry::PointCloud", "url": "class_rodin_1_1_geometry_1_1_point_cloud.html", "keys": ["rodin::geometry::pointcloud", "geometry::pointcloud", "pointcloud"]}, {"name": "Rodin::Geometry::Polytope", "url": "class_rodin_1_1_geometry_1_1_polytope.html", "keys": ["rodin::geometry::polytope", "geometry::polytope", "polytope"]}, {"name": "Rodin::Geometry::Polytope::Key", "url": "class_rodin_1_1_geometry_1_1_polytope_1_1_key.html", "keys": ["rodin::geometry::polytope::key", "geometry::polytope::key", "polytope::key", "key"]}, {"name": "Rodin::Geometry::PolytopeIterator", "url": "class_rodin_1_1_geometry_1_1_polytope_iterator.html", "keys": ["rodin::geometry::polytopeiterator", "geometry::polytopeiterator", "polytopeiterator"]}, {"name": "Rodin::Geometry::PolytopeIteratorBase", "url": "class_rodin_1_1_geometry_1_1_polytope_iterator_base.html", "keys": ["rodin::geometry::polytopeiteratorbase", "geometry::polytopeiteratorbase", "polytopeiteratorbase"]}, {"name": "Rodin::Geometry::PolytopeQuadrature", "url": "class_rodin_1_1_geometry_1_1_polytope_quadrature.html", "keys": ["rodin::geometry::polytopequadrature", "geometry::polytopequadrature", "polytopequadrature"]}, {"name": "Rodin::Geometry::PolytopeQuadratureIndex", "url": "class_rodin_1_1_geometry_1_1_polytope_quadrature_index.html", "keys": ["rodin::geometry::polytopequadratureindex", "geometry::polytopequadratureindex", "polytopequadratureindex"]}, {"name": "Rodin::Geometry::PolytopeTransformation", "url": "class_rodin_1_1_geometry_1_1_polytope_transformation.html", "keys": ["rodin::geometry::polytopetransformation", "geometry::polytopetransformation", "polytopetransformation"]}, {"name": "Rodin::Geometry::PolytopeTransformationIndex", "url": "class_rodin_1_1_geometry_1_1_polytope_transformation_index.html", "keys": ["rodin::geometry::polytopetransformationindex", "geometry::polytopetransformationindex", "polytopetransformationindex"]}, {"name": "Rodin::Geometry::SetIndexGenerator", "url": "class_rodin_1_1_geometry_1_1_set_index_generator.html", "keys": ["rodin::geometry::setindexgenerator", "geometry::setindexgenerator", "setindexgenerator"]}, {"name": "Rodin::Geometry::Shard", "url": "class_rodin_1_1_geometry_1_1_shard.html", "keys": ["rodin::geometry::shard", "geometry::shard", "shard"]}, {"name": "Rodin::Geometry::Shard::Builder", "url": "class_rodin_1_1_geometry_1_1_shard_1_1_builder.html", "keys": ["rodin::geometry::shard::builder", "geometry::shard::builder", "shard::builder", "builder"]}, {"name": "Rodin::Geometry::Sharder", "url": "class_rodin_1_1_geometry_1_1_sharder.html", "keys": ["rodin::geometry::sharder", "geometry::sharder", "sharder"]}, {"name": "Rodin::Geometry::SharderBase", "url": "class_rodin_1_1_geometry_1_1_sharder_base.html", "keys": ["rodin::geometry::sharderbase", "geometry::sharderbase", "sharderbase"]}, {"name": "Rodin::Geometry::SubMesh", "url": "class_rodin_1_1_geometry_1_1_sub_mesh.html", "keys": ["rodin::geometry::submesh", "geometry::submesh", "submesh"]}, {"name": "Rodin::Geometry::SubMeshBase", "url": "class_rodin_1_1_geometry_1_1_sub_mesh_base.html", "keys": ["rodin::geometry::submeshbase", "geometry::submeshbase", "submeshbase"]}, {"name": "Rodin::Geometry::SubMeshBuilder", "url": "class_rodin_1_1_geometry_1_1_sub_mesh_builder.html", "keys": ["rodin::geometry::submeshbuilder", "geometry::submeshbuilder", "submeshbuilder"]}, {"name": "Rodin::Geometry::VectorIndexGenerator", "url": "class_rodin_1_1_geometry_1_1_vector_index_generator.html", "keys": ["rodin::geometry::vectorindexgenerator", "geometry::vectorindexgenerator", "vectorindexgenerator"]}, {"name": "Rodin::Geometry::Vertex", "url": "class_rodin_1_1_geometry_1_1_vertex.html", "keys": ["rodin::geometry::vertex", "geometry::vertex", "vertex"]}, {"name": "Rodin::Geometry::VertexIterator", "url": "class_rodin_1_1_geometry_1_1_vertex_iterator.html", "keys": ["rodin::geometry::vertexiterator", "geometry::vertexiterator", "vertexiterator"]}, {"name": "Rodin::Heart::CCMLC2014::Numerics::DynamicSystem", "url": "class_rodin_1_1_heart_1_1_c_c_m_l_c2014_1_1_numerics_1_1_dynamic_system.html", "keys": ["rodin::heart::ccmlc2014::numerics::dynamicsystem", "heart::ccmlc2014::numerics::dynamicsystem", "ccmlc2014::numerics::dynamicsystem", "numerics::dynamicsystem", "dynamicsystem"]}, {"name": "Rodin::Heart::CCMLC2014::Physics::CavityFluxEvaluator", "url": "class_rodin_1_1_heart_1_1_c_c_m_l_c2014_1_1_physics_1_1_cavity_flux_evaluator.html", "keys": ["rodin::heart::ccmlc2014::physics::cavityfluxevaluator", "heart::ccmlc2014::physics::cavityfluxevaluator", "ccmlc2014::physics::cavityfluxevaluator", "physics::cavityfluxevaluator", "cavityfluxevaluator"]}, {"name": "Rodin::Heart::CCMLC2014::Physics::MidpointKinematics", "url": "class_rodin_1_1_heart_1_1_c_c_m_l_c2014_1_1_physics_1_1_midpoint_kinematics.html", "keys": ["rodin::heart::ccmlc2014::physics::midpointkinematics", "heart::ccmlc2014::physics::midpointkinematics", "ccmlc2014::physics::midpointkinematics", "physics::midpointkinematics", "midpointkinematics"]}, {"name": "Rodin::Heart::CCMLC2014::Physics::PassiveStressEvaluator", "url": "class_rodin_1_1_heart_1_1_c_c_m_l_c2014_1_1_physics_1_1_passive_stress_evaluator.html", "keys": ["rodin::heart::ccmlc2014::physics::passivestressevaluator", "heart::ccmlc2014::physics::passivestressevaluator", "ccmlc2014::physics::passivestressevaluator", "physics::passivestressevaluator", "passivestressevaluator"]}, {"name": "Rodin::Heart::CCMLC2014::Physics::ViscousStressEvaluator", "url": "class_rodin_1_1_heart_1_1_c_c_m_l_c2014_1_1_physics_1_1_viscous_stress_evaluator.html", "keys": ["rodin::heart::ccmlc2014::physics::viscousstressevaluator", "heart::ccmlc2014::physics::viscousstressevaluator", "ccmlc2014::physics::viscousstressevaluator", "physics::viscousstressevaluator", "viscousstressevaluator"]}, {"name": "Rodin::Heart::CCMLC2014::Physics::WindkesselOutflowEvaluator", "url": "class_rodin_1_1_heart_1_1_c_c_m_l_c2014_1_1_physics_1_1_windkessel_outflow_evaluator.html", "keys": ["rodin::heart::ccmlc2014::physics::windkesseloutflowevaluator", "heart::ccmlc2014::physics::windkesseloutflowevaluator", "ccmlc2014::physics::windkesseloutflowevaluator", "physics::windkesseloutflowevaluator", "windkesseloutflowevaluator"]}, {"name": "Rodin::Heart::CCMLC2014::Solver::StepperT", "url": "class_rodin_1_1_heart_1_1_c_c_m_l_c2014_1_1_solver_1_1_stepper_t.html", "keys": ["rodin::heart::ccmlc2014::solver::steppert", "heart::ccmlc2014::solver::steppert", "ccmlc2014::solver::steppert", "solver::steppert", "steppert"]}, {"name": "Rodin::Heart::CCMLC2014::Solver::StepperT::Problem", "url": "class_rodin_1_1_heart_1_1_c_c_m_l_c2014_1_1_solver_1_1_stepper_t_1_1_problem.html", "keys": ["rodin::heart::ccmlc2014::solver::steppert::problem", "heart::ccmlc2014::solver::steppert::problem", "ccmlc2014::solver::steppert::problem", "solver::steppert::problem", "steppert::problem", "problem"]}, {"name": "Rodin::Heart::HolzapfelReducedLaw", "url": "class_rodin_1_1_heart_1_1_holzapfel_reduced_law.html", "keys": ["rodin::heart::holzapfelreducedlaw", "heart::holzapfelreducedlaw", "holzapfelreducedlaw"]}, {"name": "Rodin::Hilbert::H1a", "url": "class_rodin_1_1_hilbert_1_1_h1a.html", "keys": ["rodin::hilbert::h1a", "hilbert::h1a", "h1a"]}, {"name": "Rodin::IO::GridFunctionLoader", "url": "class_rodin_1_1_i_o_1_1_grid_function_loader.html", "keys": ["rodin::io::gridfunctionloader", "io::gridfunctionloader", "gridfunctionloader"]}, {"name": "Rodin::IO::GridFunctionLoaderBase", "url": "class_rodin_1_1_i_o_1_1_grid_function_loader_base.html", "keys": ["rodin::io::gridfunctionloaderbase", "io::gridfunctionloaderbase", "gridfunctionloaderbase"]}, {"name": "Rodin::IO::GridFunctionPrinter", "url": "class_rodin_1_1_i_o_1_1_grid_function_printer.html", "keys": ["rodin::io::gridfunctionprinter", "io::gridfunctionprinter", "gridfunctionprinter"]}, {"name": "Rodin::IO::GridFunctionPrinterBase", "url": "class_rodin_1_1_i_o_1_1_grid_function_printer_base.html", "keys": ["rodin::io::gridfunctionprinterbase", "io::gridfunctionprinterbase", "gridfunctionprinterbase"]}, {"name": "Rodin::IO::HDF5::Handle", "url": "class_rodin_1_1_i_o_1_1_h_d_f5_1_1_handle.html", "keys": ["rodin::io::hdf5::handle", "io::hdf5::handle", "hdf5::handle", "handle"]}, {"name": "Rodin::IO::Loader", "url": "class_rodin_1_1_i_o_1_1_loader.html", "keys": ["rodin::io::loader", "io::loader", "loader"]}, {"name": "Rodin::IO::MEDIT::ParseEntity", "url": "class_rodin_1_1_i_o_1_1_m_e_d_i_t_1_1_parse_entity.html", "keys": ["rodin::io::medit::parseentity", "io::medit::parseentity", "medit::parseentity", "parseentity"]}, {"name": "Rodin::IO::MFEM::ParseMeshHeader", "url": "class_rodin_1_1_i_o_1_1_m_f_e_m_1_1_parse_mesh_header.html", "keys": ["rodin::io::mfem::parsemeshheader", "io::mfem::parsemeshheader", "mfem::parsemeshheader", "parsemeshheader"]}, {"name": "Rodin::IO::MFEM::ParseVertex", "url": "class_rodin_1_1_i_o_1_1_m_f_e_m_1_1_parse_vertex.html", "keys": ["rodin::io::mfem::parsevertex", "io::mfem::parsevertex", "mfem::parsevertex", "parsevertex"]}, {"name": "Rodin::IO::MFEM::TetrahedronNodes", "url": "class_rodin_1_1_i_o_1_1_m_f_e_m_1_1_tetrahedron_nodes.html", "keys": ["rodin::io::mfem::tetrahedronnodes", "io::mfem::tetrahedronnodes", "mfem::tetrahedronnodes", "tetrahedronnodes"]}, {"name": "Rodin::IO::MFEM::TriangleNodes", "url": "class_rodin_1_1_i_o_1_1_m_f_e_m_1_1_triangle_nodes.html", "keys": ["rodin::io::mfem::trianglenodes", "io::mfem::trianglenodes", "mfem::trianglenodes", "trianglenodes"]}, {"name": "Rodin::IO::MFEM::VandermondeTetrahedron", "url": "class_rodin_1_1_i_o_1_1_m_f_e_m_1_1_vandermonde_tetrahedron.html", "keys": ["rodin::io::mfem::vandermondetetrahedron", "io::mfem::vandermondetetrahedron", "mfem::vandermondetetrahedron", "vandermondetetrahedron"]}, {"name": "Rodin::IO::MFEM::VandermondeTriangle", "url": "class_rodin_1_1_i_o_1_1_m_f_e_m_1_1_vandermonde_triangle.html", "keys": ["rodin::io::mfem::vandermondetriangle", "io::mfem::vandermondetriangle", "mfem::vandermondetriangle", "vandermondetriangle"]}, {"name": "Rodin::IO::MFEM::WedgeChange", "url": "class_rodin_1_1_i_o_1_1_m_f_e_m_1_1_wedge_change.html", "keys": ["rodin::io::mfem::wedgechange", "io::mfem::wedgechange", "mfem::wedgechange", "wedgechange"]}, {"name": "Rodin::IO::MFEM::WedgeNodes", "url": "class_rodin_1_1_i_o_1_1_m_f_e_m_1_1_wedge_nodes.html", "keys": ["rodin::io::mfem::wedgenodes", "io::mfem::wedgenodes", "mfem::wedgenodes", "wedgenodes"]}, {"name": "Rodin::IO::MeshLoader", "url": "class_rodin_1_1_i_o_1_1_mesh_loader.html", "keys": ["rodin::io::meshloader", "io::meshloader", "meshloader"]}, {"name": "Rodin::IO::MeshLoaderBase", "url": "class_rodin_1_1_i_o_1_1_mesh_loader_base.html", "keys": ["rodin::io::meshloaderbase", "io::meshloaderbase", "meshloaderbase"]}, {"name": "Rodin::IO::MeshPrinter", "url": "class_rodin_1_1_i_o_1_1_mesh_printer.html", "keys": ["rodin::io::meshprinter", "io::meshprinter", "meshprinter"]}, {"name": "Rodin::IO::MeshPrinterBase", "url": "class_rodin_1_1_i_o_1_1_mesh_printer_base.html", "keys": ["rodin::io::meshprinterbase", "io::meshprinterbase", "meshprinterbase"]}, {"name": "Rodin::IO::Printer", "url": "class_rodin_1_1_i_o_1_1_printer.html", "keys": ["rodin::io::printer", "io::printer", "printer"]}, {"name": "Rodin::IO::XDMF", "url": "class_rodin_1_1_i_o_1_1_x_d_m_f.html", "keys": ["rodin::io::xdmf", "io::xdmf", "xdmf"]}, {"name": "Rodin::IO::XDMF::Grid", "url": "class_rodin_1_1_i_o_1_1_x_d_m_f_1_1_grid.html", "keys": ["rodin::io::xdmf::grid", "io::xdmf::grid", "xdmf::grid", "grid"]}, {"name": "Rodin::Identifiable", "url": "class_rodin_1_1_identifiable.html", "keys": ["rodin::identifiable", "identifiable"]}, {"name": "Rodin::Location::AABB", "url": "class_rodin_1_1_location_1_1_a_a_b_b.html", "keys": ["rodin::location::aabb", "location::aabb", "aabb"]}, {"name": "Rodin::MMG::Adapt", "url": "class_rodin_1_1_m_m_g_1_1_adapt.html", "keys": ["rodin::mmg::adapt", "mmg::adapt", "adapt"]}, {"name": "Rodin::MMG::Base", "url": "class_rodin_1_1_m_m_g_1_1_base.html", "keys": ["rodin::mmg::base", "mmg::base", "base"]}, {"name": "Rodin::MMG::LevelSetDiscretizer", "url": "class_rodin_1_1_m_m_g_1_1_level_set_discretizer.html", "keys": ["rodin::mmg::levelsetdiscretizer", "mmg::levelsetdiscretizer", "levelsetdiscretizer"]}, {"name": "Rodin::MMG::MMGException", "url": "class_rodin_1_1_m_m_g_1_1_m_m_g_exception.html", "keys": ["rodin::mmg::mmgexception", "mmg::mmgexception", "mmgexception"]}, {"name": "Rodin::MMG::Mesh", "url": "class_rodin_1_1_m_m_g_1_1_mesh.html", "keys": ["rodin::mmg::mesh", "mmg::mesh", "mesh"]}, {"name": "Rodin::MMG::Mesh::Builder", "url": "class_rodin_1_1_m_m_g_1_1_mesh_1_1_builder.html", "keys": ["rodin::mmg::mesh::builder", "mmg::mesh::builder", "mesh::builder", "builder"]}, {"name": "Rodin::MMG::MeshLoader", "url": "class_rodin_1_1_m_m_g_1_1_mesh_loader.html", "keys": ["rodin::mmg::meshloader", "mmg::meshloader", "meshloader"]}, {"name": "Rodin::MMG::MeshPrinter", "url": "class_rodin_1_1_m_m_g_1_1_mesh_printer.html", "keys": ["rodin::mmg::meshprinter", "mmg::meshprinter", "meshprinter"]}, {"name": "Rodin::MMG::Optimizer", "url": "class_rodin_1_1_m_m_g_1_1_optimizer.html", "keys": ["rodin::mmg::optimizer", "mmg::optimizer", "optimizer"]}, {"name": "Rodin::Math::Deg", "url": "class_rodin_1_1_math_1_1_deg.html", "keys": ["rodin::math::deg", "math::deg", "deg"]}, {"name": "Rodin::Math::FixedSizeMatrix", "url": "class_rodin_1_1_math_1_1_fixed_size_matrix.html", "keys": ["rodin::math::fixedsizematrix", "math::fixedsizematrix", "fixedsizematrix"]}, {"name": "Rodin::Math::FixedSizeVector", "url": "class_rodin_1_1_math_1_1_fixed_size_vector.html", "keys": ["rodin::math::fixedsizevector", "math::fixedsizevector", "fixedsizevector"]}, {"name": "Rodin::Math::LinearSystem", "url": "class_rodin_1_1_math_1_1_linear_system.html", "keys": ["rodin::math::linearsystem", "math::linearsystem", "linearsystem"]}, {"name": "Rodin::Math::LinearSystemBase", "url": "class_rodin_1_1_math_1_1_linear_system_base.html", "keys": ["rodin::math::linearsystembase", "math::linearsystembase", "linearsystembase"]}, {"name": "Rodin::Math::Matrix", "url": "class_rodin_1_1_math_1_1_matrix.html", "keys": ["rodin::math::matrix", "math::matrix", "matrix"]}, {"name": "Rodin::Math::Rad", "url": "class_rodin_1_1_math_1_1_rad.html", "keys": ["rodin::math::rad", "math::rad", "rad"]}, {"name": "Rodin::Math::RootFinding::NewtonRaphson", "url": "class_rodin_1_1_math_1_1_root_finding_1_1_newton_raphson.html", "keys": ["rodin::math::rootfinding::newtonraphson", "math::rootfinding::newtonraphson", "rootfinding::newtonraphson", "newtonraphson"]}, {"name": "Rodin::Math::SpatialMatrix", "url": "class_rodin_1_1_math_1_1_spatial_matrix.html", "keys": ["rodin::math::spatialmatrix", "math::spatialmatrix", "spatialmatrix"]}, {"name": "Rodin::Math::SpatialTensor", "url": "class_rodin_1_1_math_1_1_spatial_tensor.html", "keys": ["rodin::math::spatialtensor", "math::spatialtensor", "spatialtensor"]}, {"name": "Rodin::Math::SpatialVector", "url": "class_rodin_1_1_math_1_1_spatial_vector.html", "keys": ["rodin::math::spatialvector", "math::spatialvector", "spatialvector"]}, {"name": "Rodin::Math::Unit", "url": "class_rodin_1_1_math_1_1_unit.html", "keys": ["rodin::math::unit", "math::unit", "unit"]}, {"name": "Rodin::Math::Vector", "url": "class_rodin_1_1_math_1_1_vector.html", "keys": ["rodin::math::vector", "math::vector", "vector"]}, {"name": "Rodin::Moveable", "url": "class_rodin_1_1_moveable.html", "keys": ["rodin::moveable", "moveable"]}, {"name": "Rodin::PETSc::Assembly::ConstraintMap", "url": "class_rodin_1_1_p_e_t_sc_1_1_assembly_1_1_constraint_map.html", "keys": ["rodin::petsc::assembly::constraintmap", "petsc::assembly::constraintmap", "assembly::constraintmap", "constraintmap"]}, {"name": "Rodin::PETSc::Assembly::MatrixSetup", "url": "class_rodin_1_1_p_e_t_sc_1_1_assembly_1_1_matrix_setup.html", "keys": ["rodin::petsc::assembly::matrixsetup", "petsc::assembly::matrixsetup", "assembly::matrixsetup", "matrixsetup"]}, {"name": "Rodin::PETSc::Assembly::VectorSetup", "url": "class_rodin_1_1_p_e_t_sc_1_1_assembly_1_1_vector_setup.html", "keys": ["rodin::petsc::assembly::vectorsetup", "petsc::assembly::vectorsetup", "assembly::vectorsetup", "vectorsetup"]}, {"name": "Rodin::PETSc::Object", "url": "class_rodin_1_1_p_e_t_sc_1_1_object.html", "keys": ["rodin::petsc::object", "petsc::object", "object"]}, {"name": "Rodin::PETSc::Variational::GridFunction", "url": "class_rodin_1_1_p_e_t_sc_1_1_variational_1_1_grid_function.html", "keys": ["rodin::petsc::variational::gridfunction", "petsc::variational::gridfunction", "variational::gridfunction", "gridfunction"]}, {"name": "Rodin::PETSc::Variational::TestFunction", "url": "class_rodin_1_1_p_e_t_sc_1_1_variational_1_1_test_function.html", "keys": ["rodin::petsc::variational::testfunction", "petsc::variational::testfunction", "variational::testfunction", "testfunction"]}, {"name": "Rodin::PETSc::Variational::TrialFunction", "url": "class_rodin_1_1_p_e_t_sc_1_1_variational_1_1_trial_function.html", "keys": ["rodin::petsc::variational::trialfunction", "petsc::variational::trialfunction", "variational::trialfunction", "trialfunction"]}, {"name": "Rodin::Pair", "url": "class_rodin_1_1_pair.html", "keys": ["rodin::pair", "pair"]}, {"name": "Rodin::QF::Centroid", "url": "class_rodin_1_1_q_f_1_1_centroid.html", "keys": ["rodin::qf::centroid", "qf::centroid", "centroid"]}, {"name": "Rodin::QF::GaussLegendre", "url": "class_rodin_1_1_q_f_1_1_gauss_legendre.html", "keys": ["rodin::qf::gausslegendre", "qf::gausslegendre", "gausslegendre"]}, {"name": "Rodin::QF::GaussLobatto", "url": "class_rodin_1_1_q_f_1_1_gauss_lobatto.html", "keys": ["rodin::qf::gausslobatto", "qf::gausslobatto", "gausslobatto"]}, {"name": "Rodin::QF::GrundmannMoller", "url": "class_rodin_1_1_q_f_1_1_grundmann_moller.html", "keys": ["rodin::qf::grundmannmoller", "qf::grundmannmoller", "grundmannmoller"]}, {"name": "Rodin::QF::QuadratureFormula", "url": "class_rodin_1_1_q_f_1_1_quadrature_formula.html", "keys": ["rodin::qf::quadratureformula", "qf::quadratureformula", "quadratureformula"]}, {"name": "Rodin::QF::QuadratureFormulaBase", "url": "class_rodin_1_1_q_f_1_1_quadrature_formula_base.html", "keys": ["rodin::qf::quadratureformulabase", "qf::quadratureformulabase", "quadratureformulabase"]}, {"name": "Rodin::QF::TensorProduct", "url": "class_rodin_1_1_q_f_1_1_tensor_product.html", "keys": ["rodin::qf::tensorproduct", "qf::tensorproduct", "tensorproduct"]}, {"name": "Rodin::QF::WitherdenVincent", "url": "class_rodin_1_1_q_f_1_1_witherden_vincent.html", "keys": ["rodin::qf::witherdenvincent", "qf::witherdenvincent", "witherdenvincent"]}, {"name": "Rodin::QF::XiaoGimbutas", "url": "class_rodin_1_1_q_f_1_1_xiao_gimbutas.html", "keys": ["rodin::qf::xiaogimbutas", "qf::xiaogimbutas", "xiaogimbutas"]}, {"name": "Rodin::Scotch::Exception", "url": "class_rodin_1_1_scotch_1_1_exception.html", "keys": ["rodin::scotch::exception", "scotch::exception", "exception"]}, {"name": "Rodin::Scotch::Partitioner", "url": "class_rodin_1_1_scotch_1_1_partitioner.html", "keys": ["rodin::scotch::partitioner", "scotch::partitioner", "partitioner"]}, {"name": "Rodin::Solid::ActiveContraction", "url": "class_rodin_1_1_solid_1_1_active_contraction.html", "keys": ["rodin::solid::activecontraction", "solid::activecontraction", "activecontraction"]}, {"name": "Rodin::Solid::ActiveFiberLaw", "url": "class_rodin_1_1_solid_1_1_active_fiber_law.html", "keys": ["rodin::solid::activefiberlaw", "solid::activefiberlaw", "activefiberlaw"]}, {"name": "Rodin::Solid::CauchyStress", "url": "class_rodin_1_1_solid_1_1_cauchy_stress.html", "keys": ["rodin::solid::cauchystress", "solid::cauchystress", "cauchystress"]}, {"name": "Rodin::Solid::ConstitutivePoint", "url": "class_rodin_1_1_solid_1_1_constitutive_point.html", "keys": ["rodin::solid::constitutivepoint", "solid::constitutivepoint", "constitutivepoint"]}, {"name": "Rodin::Solid::FiberInvariants", "url": "class_rodin_1_1_solid_1_1_fiber_invariants.html", "keys": ["rodin::solid::fiberinvariants", "solid::fiberinvariants", "fiberinvariants"]}, {"name": "Rodin::Solid::FiberKinematics", "url": "class_rodin_1_1_solid_1_1_fiber_kinematics.html", "keys": ["rodin::solid::fiberkinematics", "solid::fiberkinematics", "fiberkinematics"]}, {"name": "Rodin::Solid::FirstPiolaKirchhoffStress", "url": "class_rodin_1_1_solid_1_1_first_piola_kirchhoff_stress.html", "keys": ["rodin::solid::firstpiolakirchhoffstress", "solid::firstpiolakirchhoffstress", "firstpiolakirchhoffstress"]}, {"name": "Rodin::Solid::FollowerPressureForce", "url": "class_rodin_1_1_solid_1_1_follower_pressure_force.html", "keys": ["rodin::solid::followerpressureforce", "solid::followerpressureforce", "followerpressureforce"]}, {"name": "Rodin::Solid::FollowerPressureTangent", "url": "class_rodin_1_1_solid_1_1_follower_pressure_tangent.html", "keys": ["rodin::solid::followerpressuretangent", "solid::followerpressuretangent", "followerpressuretangent"]}, {"name": "Rodin::Solid::GreenLagrangeStrain", "url": "class_rodin_1_1_solid_1_1_green_lagrange_strain.html", "keys": ["rodin::solid::greenlagrangestrain", "solid::greenlagrangestrain", "greenlagrangestrain"]}, {"name": "Rodin::Solid::HolzapfelOgden", "url": "class_rodin_1_1_solid_1_1_holzapfel_ogden.html", "keys": ["rodin::solid::holzapfelogden", "solid::holzapfelogden", "holzapfelogden"]}, {"name": "Rodin::Solid::Hooke", "url": "class_rodin_1_1_solid_1_1_hooke.html", "keys": ["rodin::solid::hooke", "solid::hooke", "hooke"]}, {"name": "Rodin::Solid::HyperElasticLaw", "url": "class_rodin_1_1_solid_1_1_hyper_elastic_law.html", "keys": ["rodin::solid::hyperelasticlaw", "solid::hyperelasticlaw", "hyperelasticlaw"]}, {"name": "Rodin::Solid::Input", "url": "class_rodin_1_1_solid_1_1_input.html", "keys": ["rodin::solid::input", "solid::input", "input"]}, {"name": "Rodin::Solid::InternalVirtualWork", "url": "class_rodin_1_1_solid_1_1_internal_virtual_work.html", "keys": ["rodin::solid::internalvirtualwork", "solid::internalvirtualwork", "internalvirtualwork"]}, {"name": "Rodin::Solid::InternalVirtualWorkResidual", "url": "class_rodin_1_1_solid_1_1_internal_virtual_work_residual.html", "keys": ["rodin::solid::internalvirtualworkresidual", "solid::internalvirtualworkresidual", "internalvirtualworkresidual"]}, {"name": "Rodin::Solid::InternalVirtualWorkResidualP", "url": "class_rodin_1_1_solid_1_1_internal_virtual_work_residual_p.html", "keys": ["rodin::solid::internalvirtualworkresidualp", "solid::internalvirtualworkresidualp", "internalvirtualworkresidualp"]}, {"name": "Rodin::Solid::InternalVirtualWorkTangent", "url": "class_rodin_1_1_solid_1_1_internal_virtual_work_tangent.html", "keys": ["rodin::solid::internalvirtualworktangent", "solid::internalvirtualworktangent", "internalvirtualworktangent"]}, {"name": "Rodin::Solid::InternalVirtualWorkTangentPU", "url": "class_rodin_1_1_solid_1_1_internal_virtual_work_tangent_p_u.html", "keys": ["rodin::solid::internalvirtualworktangentpu", "solid::internalvirtualworktangentpu", "internalvirtualworktangentpu"]}, {"name": "Rodin::Solid::InternalVirtualWorkTangentUP", "url": "class_rodin_1_1_solid_1_1_internal_virtual_work_tangent_u_p.html", "keys": ["rodin::solid::internalvirtualworktangentup", "solid::internalvirtualworktangentup", "internalvirtualworktangentup"]}, {"name": "Rodin::Solid::IsotropicInvariants", "url": "class_rodin_1_1_solid_1_1_isotropic_invariants.html", "keys": ["rodin::solid::isotropicinvariants", "solid::isotropicinvariants", "isotropicinvariants"]}, {"name": "Rodin::Solid::KinematicState", "url": "class_rodin_1_1_solid_1_1_kinematic_state.html", "keys": ["rodin::solid::kinematicstate", "solid::kinematicstate", "kinematicstate"]}, {"name": "Rodin::Solid::MooneyRivlin", "url": "class_rodin_1_1_solid_1_1_mooney_rivlin.html", "keys": ["rodin::solid::mooneyrivlin", "solid::mooneyrivlin", "mooneyrivlin"]}, {"name": "Rodin::Solid::NeoHookean", "url": "class_rodin_1_1_solid_1_1_neo_hookean.html", "keys": ["rodin::solid::neohookean", "solid::neohookean", "neohookean"]}, {"name": "Rodin::Solid::Output", "url": "class_rodin_1_1_solid_1_1_output.html", "keys": ["rodin::solid::output", "solid::output", "output"]}, {"name": "Rodin::Solid::SaintVenantKirchhoff", "url": "class_rodin_1_1_solid_1_1_saint_venant_kirchhoff.html", "keys": ["rodin::solid::saintvenantkirchhoff", "solid::saintvenantkirchhoff", "saintvenantkirchhoff"]}, {"name": "Rodin::Solid::Yeoh", "url": "class_rodin_1_1_solid_1_1_yeoh.html", "keys": ["rodin::solid::yeoh", "solid::yeoh", "yeoh"]}, {"name": "Rodin::Solver::BiCGSTAB", "url": "class_rodin_1_1_solver_1_1_bi_c_g_s_t_a_b.html", "keys": ["rodin::solver::bicgstab", "solver::bicgstab", "bicgstab"]}, {"name": "Rodin::Solver::CG", "url": "class_rodin_1_1_solver_1_1_c_g.html", "keys": ["rodin::solver::cg", "solver::cg", "cg"]}, {"name": "Rodin::Solver::DGMRES", "url": "class_rodin_1_1_solver_1_1_d_g_m_r_e_s.html", "keys": ["rodin::solver::dgmres", "solver::dgmres", "dgmres"]}, {"name": "Rodin::Solver::EigenSolver", "url": "class_rodin_1_1_solver_1_1_eigen_solver.html", "keys": ["rodin::solver::eigensolver", "solver::eigensolver", "eigensolver"]}, {"name": "Rodin::Solver::GMRES", "url": "class_rodin_1_1_solver_1_1_g_m_r_e_s.html", "keys": ["rodin::solver::gmres", "solver::gmres", "gmres"]}, {"name": "Rodin::Solver::HouseholderQR", "url": "class_rodin_1_1_solver_1_1_householder_q_r.html", "keys": ["rodin::solver::householderqr", "solver::householderqr", "householderqr"]}, {"name": "Rodin::Solver::IDRS", "url": "class_rodin_1_1_solver_1_1_i_d_r_s.html", "keys": ["rodin::solver::idrs", "solver::idrs", "idrs"]}, {"name": "Rodin::Solver::IDRSTABL", "url": "class_rodin_1_1_solver_1_1_i_d_r_s_t_a_b_l.html", "keys": ["rodin::solver::idrstabl", "solver::idrstabl", "idrstabl"]}, {"name": "Rodin::Solver::KSP", "url": "class_rodin_1_1_solver_1_1_k_s_p.html", "keys": ["rodin::solver::ksp", "solver::ksp", "ksp"]}, {"name": "Rodin::Solver::LDLT", "url": "class_rodin_1_1_solver_1_1_l_d_l_t.html", "keys": ["rodin::solver::ldlt", "solver::ldlt", "ldlt"]}, {"name": "Rodin::Solver::LeastSquaresCG", "url": "class_rodin_1_1_solver_1_1_least_squares_c_g.html", "keys": ["rodin::solver::leastsquarescg", "solver::leastsquarescg", "leastsquarescg"]}, {"name": "Rodin::Solver::LinearSolverBase", "url": "class_rodin_1_1_solver_1_1_linear_solver_base.html", "keys": ["rodin::solver::linearsolverbase", "solver::linearsolverbase", "linearsolverbase"]}, {"name": "Rodin::Solver::MINRES", "url": "class_rodin_1_1_solver_1_1_m_i_n_r_e_s.html", "keys": ["rodin::solver::minres", "solver::minres", "minres"]}, {"name": "Rodin::Solver::NewtonSolver", "url": "class_rodin_1_1_solver_1_1_newton_solver.html", "keys": ["rodin::solver::newtonsolver", "solver::newtonsolver", "newtonsolver"]}, {"name": "Rodin::Solver::NewtonSolverBase", "url": "class_rodin_1_1_solver_1_1_newton_solver_base.html", "keys": ["rodin::solver::newtonsolverbase", "solver::newtonsolverbase", "newtonsolverbase"]}, {"name": "Rodin::Solver::PartialPivLU", "url": "class_rodin_1_1_solver_1_1_partial_piv_l_u.html", "keys": ["rodin::solver::partialpivlu", "solver::partialpivlu", "partialpivlu"]}, {"name": "Rodin::Solver::SNES", "url": "class_rodin_1_1_solver_1_1_s_n_e_s.html", "keys": ["rodin::solver::snes", "solver::snes", "snes"]}, {"name": "Rodin::Solver::SimplicialLDLT", "url": "class_rodin_1_1_solver_1_1_simplicial_l_d_l_t.html", "keys": ["rodin::solver::simplicialldlt", "solver::simplicialldlt", "simplicialldlt"]}, {"name": "Rodin::Solver::SimplicialLLT", "url": "class_rodin_1_1_solver_1_1_simplicial_l_l_t.html", "keys": ["rodin::solver::simplicialllt", "solver::simplicialllt", "simplicialllt"]}, {"name": "Rodin::Solver::SparseLU", "url": "class_rodin_1_1_solver_1_1_sparse_l_u.html", "keys": ["rodin::solver::sparselu", "solver::sparselu", "sparselu"]}, {"name": "Rodin::Solver::SparseQR", "url": "class_rodin_1_1_solver_1_1_sparse_q_r.html", "keys": ["rodin::solver::sparseqr", "solver::sparseqr", "sparseqr"]}, {"name": "Rodin::Test::FDProbe", "url": "class_rodin_1_1_test_1_1_f_d_probe.html", "keys": ["rodin::test::fdprobe", "test::fdprobe", "fdprobe"]}, {"name": "Rodin::Test::Random::PointOnTriangle", "url": "class_rodin_1_1_test_1_1_random_1_1_point_on_triangle.html", "keys": ["rodin::test::random::pointontriangle", "test::random::pointontriangle", "random::pointontriangle", "pointontriangle"]}, {"name": "Rodin::Test::Random::RandomFloat", "url": "class_rodin_1_1_test_1_1_random_1_1_random_float.html", "keys": ["rodin::test::random::randomfloat", "test::random::randomfloat", "random::randomfloat", "randomfloat"]}, {"name": "Rodin::Test::Random::RandomInteger", "url": "class_rodin_1_1_test_1_1_random_1_1_random_integer.html", "keys": ["rodin::test::random::randominteger", "test::random::randominteger", "random::randominteger", "randominteger"]}, {"name": "Rodin::Test::Utility::ExpectedResultTable", "url": "class_rodin_1_1_test_1_1_utility_1_1_expected_result_table.html", "keys": ["rodin::test::utility::expectedresulttable", "test::utility::expectedresulttable", "utility::expectedresulttable", "expectedresulttable"]}, {"name": "Rodin::Test::Utility::ExpectedResultTable::ExpectedResult", "url": "class_rodin_1_1_test_1_1_utility_1_1_expected_result_table_1_1_expected_result.html", "keys": ["rodin::test::utility::expectedresulttable::expectedresult", "test::utility::expectedresulttable::expectedresult", "utility::expectedresulttable::expectedresult", "expectedresulttable::expectedresult", "expectedresult"]}, {"name": "Rodin::Tuple", "url": "class_rodin_1_1_tuple.html", "keys": ["rodin::tuple", "tuple"]}, {"name": "Rodin::Utility::Cast", "url": "class_rodin_1_1_utility_1_1_cast.html", "keys": ["rodin::utility::cast", "utility::cast", "cast"]}, {"name": "Rodin::Utility::IntegerSequence", "url": "class_rodin_1_1_utility_1_1_integer_sequence.html", "keys": ["rodin::utility::integersequence", "utility::integersequence", "integersequence"]}, {"name": "Rodin::Utility::OptionalReference", "url": "class_rodin_1_1_utility_1_1_optional_reference.html", "keys": ["rodin::utility::optionalreference", "utility::optionalreference", "optionalreference"]}, {"name": "Rodin::Utility::ParameterPack", "url": "class_rodin_1_1_utility_1_1_parameter_pack.html", "keys": ["rodin::utility::parameterpack", "utility::parameterpack", "parameterpack"]}, {"name": "Rodin::Utility::Wrap", "url": "class_rodin_1_1_utility_1_1_wrap.html", "keys": ["rodin::utility::wrap", "utility::wrap", "wrap"]}, {"name": "Rodin::Variational::AND", "url": "class_rodin_1_1_variational_1_1_a_n_d.html", "keys": ["rodin::variational::and", "variational::and", "and"]}, {"name": "Rodin::Variational::Abs", "url": "class_rodin_1_1_variational_1_1_abs.html", "keys": ["rodin::variational::abs", "variational::abs", "abs"]}, {"name": "Rodin::Variational::Average", "url": "class_rodin_1_1_variational_1_1_average.html", "keys": ["rodin::variational::average", "variational::average", "average"]}, {"name": "Rodin::Variational::BernsteinPyramid", "url": "class_rodin_1_1_variational_1_1_bernstein_pyramid.html", "keys": ["rodin::variational::bernsteinpyramid", "variational::bernsteinpyramid", "bernsteinpyramid"]}, {"name": "Rodin::Variational::BilinearForm", "url": "class_rodin_1_1_variational_1_1_bilinear_form.html", "keys": ["rodin::variational::bilinearform", "variational::bilinearform", "bilinearform"]}, {"name": "Rodin::Variational::BilinearFormBase", "url": "class_rodin_1_1_variational_1_1_bilinear_form_base.html", "keys": ["rodin::variational::bilinearformbase", "variational::bilinearformbase", "bilinearformbase"]}, {"name": "Rodin::Variational::BilinearFormIntegratorBase", "url": "class_rodin_1_1_variational_1_1_bilinear_form_integrator_base.html", "keys": ["rodin::variational::bilinearformintegratorbase", "variational::bilinearformintegratorbase", "bilinearformintegratorbase"]}, {"name": "Rodin::Variational::BooleanFunction", "url": "class_rodin_1_1_variational_1_1_boolean_function.html", "keys": ["rodin::variational::booleanfunction", "variational::booleanfunction", "booleanfunction"]}, {"name": "Rodin::Variational::BooleanFunctionBase", "url": "class_rodin_1_1_variational_1_1_boolean_function_base.html", "keys": ["rodin::variational::booleanfunctionbase", "variational::booleanfunctionbase", "booleanfunctionbase"]}, {"name": "Rodin::Variational::BoundaryIntegral", "url": "class_rodin_1_1_variational_1_1_boundary_integral.html", "keys": ["rodin::variational::boundaryintegral", "variational::boundaryintegral", "boundaryintegral"]}, {"name": "Rodin::Variational::BoundaryNormal", "url": "class_rodin_1_1_variational_1_1_boundary_normal.html", "keys": ["rodin::variational::boundarynormal", "variational::boundarynormal", "boundarynormal"]}, {"name": "Rodin::Variational::ComplexFunction", "url": "class_rodin_1_1_variational_1_1_complex_function.html", "keys": ["rodin::variational::complexfunction", "variational::complexfunction", "complexfunction"]}, {"name": "Rodin::Variational::ComplexFunctionBase", "url": "class_rodin_1_1_variational_1_1_complex_function_base.html", "keys": ["rodin::variational::complexfunctionbase", "variational::complexfunctionbase", "complexfunctionbase"]}, {"name": "Rodin::Variational::Component", "url": "class_rodin_1_1_variational_1_1_component.html", "keys": ["rodin::variational::component", "variational::component", "component"]}, {"name": "Rodin::Variational::Composition", "url": "class_rodin_1_1_variational_1_1_composition.html", "keys": ["rodin::variational::composition", "variational::composition", "composition"]}, {"name": "Rodin::Variational::Conjugate", "url": "class_rodin_1_1_variational_1_1_conjugate.html", "keys": ["rodin::variational::conjugate", "variational::conjugate", "conjugate"]}, {"name": "Rodin::Variational::Cos", "url": "class_rodin_1_1_variational_1_1_cos.html", "keys": ["rodin::variational::cos", "variational::cos", "cos"]}, {"name": "Rodin::Variational::Cosh", "url": "class_rodin_1_1_variational_1_1_cosh.html", "keys": ["rodin::variational::cosh", "variational::cosh", "cosh"]}, {"name": "Rodin::Variational::DOFMapping", "url": "class_rodin_1_1_variational_1_1_d_o_f_mapping.html", "keys": ["rodin::variational::dofmapping", "variational::dofmapping", "dofmapping"]}, {"name": "Rodin::Variational::DefaultBoundaryPolicy", "url": "class_rodin_1_1_variational_1_1_default_boundary_policy.html", "keys": ["rodin::variational::defaultboundarypolicy", "variational::defaultboundarypolicy", "defaultboundarypolicy"]}, {"name": "Rodin::Variational::DenseProblem", "url": "class_rodin_1_1_variational_1_1_dense_problem.html", "keys": ["rodin::variational::denseproblem", "variational::denseproblem", "denseproblem"]}, {"name": "Rodin::Variational::Derivative", "url": "class_rodin_1_1_variational_1_1_derivative.html", "keys": ["rodin::variational::derivative", "variational::derivative", "derivative"]}, {"name": "Rodin::Variational::DerivativeBase", "url": "class_rodin_1_1_variational_1_1_derivative_base.html", "keys": ["rodin::variational::derivativebase", "variational::derivativebase", "derivativebase"]}, {"name": "Rodin::Variational::DiffusionForm", "url": "class_rodin_1_1_variational_1_1_diffusion_form.html", "keys": ["rodin::variational::diffusionform", "variational::diffusionform", "diffusionform"]}, {"name": "Rodin::Variational::DiffusionForm::Kernel", "url": "class_rodin_1_1_variational_1_1_diffusion_form_1_1_kernel.html", "keys": ["rodin::variational::diffusionform::kernel", "variational::diffusionform::kernel", "diffusionform::kernel", "kernel"]}, {"name": "Rodin::Variational::DirichletBC", "url": "class_rodin_1_1_variational_1_1_dirichlet_b_c.html", "keys": ["rodin::variational::dirichletbc", "variational::dirichletbc", "dirichletbc"]}, {"name": "Rodin::Variational::DirichletBCBase", "url": "class_rodin_1_1_variational_1_1_dirichlet_b_c_base.html", "keys": ["rodin::variational::dirichletbcbase", "variational::dirichletbcbase", "dirichletbcbase"]}, {"name": "Rodin::Variational::Div", "url": "class_rodin_1_1_variational_1_1_div.html", "keys": ["rodin::variational::div", "variational::div", "div"]}, {"name": "Rodin::Variational::DivBase", "url": "class_rodin_1_1_variational_1_1_div_base.html", "keys": ["rodin::variational::divbase", "variational::divbase", "divbase"]}, {"name": "Rodin::Variational::Division", "url": "class_rodin_1_1_variational_1_1_division.html", "keys": ["rodin::variational::division", "variational::division", "division"]}, {"name": "Rodin::Variational::Dot", "url": "class_rodin_1_1_variational_1_1_dot.html", "keys": ["rodin::variational::dot", "variational::dot", "dot"]}, {"name": "Rodin::Variational::DubinerTetrahedron", "url": "class_rodin_1_1_variational_1_1_dubiner_tetrahedron.html", "keys": ["rodin::variational::dubinertetrahedron", "variational::dubinertetrahedron", "dubinertetrahedron"]}, {"name": "Rodin::Variational::DubinerTriangle", "url": "class_rodin_1_1_variational_1_1_dubiner_triangle.html", "keys": ["rodin::variational::dubinertriangle", "variational::dubinertriangle", "dubinertriangle"]}, {"name": "Rodin::Variational::EQ", "url": "class_rodin_1_1_variational_1_1_e_q.html", "keys": ["rodin::variational::eq", "variational::eq", "eq"]}, {"name": "Rodin::Variational::Exp", "url": "class_rodin_1_1_variational_1_1_exp.html", "keys": ["rodin::variational::exp", "variational::exp", "exp"]}, {"name": "Rodin::Variational::F::X", "url": "class_rodin_1_1_variational_1_1_f_1_1_x.html", "keys": ["rodin::variational::f::x", "variational::f::x", "f::x", "x"]}, {"name": "Rodin::Variational::F::Y", "url": "class_rodin_1_1_variational_1_1_f_1_1_y.html", "keys": ["rodin::variational::f::y", "variational::f::y", "f::y", "y"]}, {"name": "Rodin::Variational::F::Z", "url": "class_rodin_1_1_variational_1_1_f_1_1_z.html", "keys": ["rodin::variational::f::z", "variational::f::z", "f::z", "z"]}, {"name": "Rodin::Variational::FaceIntegral", "url": "class_rodin_1_1_variational_1_1_face_integral.html", "keys": ["rodin::variational::faceintegral", "variational::faceintegral", "faceintegral"]}, {"name": "Rodin::Variational::FaceNormal", "url": "class_rodin_1_1_variational_1_1_face_normal.html", "keys": ["rodin::variational::facenormal", "variational::facenormal", "facenormal"]}, {"name": "Rodin::Variational::FeketeTetrahedron", "url": "class_rodin_1_1_variational_1_1_fekete_tetrahedron.html", "keys": ["rodin::variational::feketetetrahedron", "variational::feketetetrahedron", "feketetetrahedron"]}, {"name": "Rodin::Variational::FeketeTriangle", "url": "class_rodin_1_1_variational_1_1_fekete_triangle.html", "keys": ["rodin::variational::feketetriangle", "variational::feketetriangle", "feketetriangle"]}, {"name": "Rodin::Variational::FiniteElementBase", "url": "class_rodin_1_1_variational_1_1_finite_element_base.html", "keys": ["rodin::variational::finiteelementbase", "variational::finiteelementbase", "finiteelementbase"]}, {"name": "Rodin::Variational::FiniteElementSpace", "url": "class_rodin_1_1_variational_1_1_finite_element_space.html", "keys": ["rodin::variational::finiteelementspace", "variational::finiteelementspace", "finiteelementspace"]}, {"name": "Rodin::Variational::FiniteElementSpaceBase", "url": "class_rodin_1_1_variational_1_1_finite_element_space_base.html", "keys": ["rodin::variational::finiteelementspacebase", "variational::finiteelementspacebase", "finiteelementspacebase"]}, {"name": "Rodin::Variational::FiniteElementSpacePullbackBase", "url": "class_rodin_1_1_variational_1_1_finite_element_space_pullback_base.html", "keys": ["rodin::variational::finiteelementspacepullbackbase", "variational::finiteelementspacepullbackbase", "finiteelementspacepullbackbase"]}, {"name": "Rodin::Variational::FiniteElementSpacePushforwardBase", "url": "class_rodin_1_1_variational_1_1_finite_element_space_pushforward_base.html", "keys": ["rodin::variational::finiteelementspacepushforwardbase", "variational::finiteelementspacepushforwardbase", "finiteelementspacepushforwardbase"]}, {"name": "Rodin::Variational::Flow", "url": "class_rodin_1_1_variational_1_1_flow.html", "keys": ["rodin::variational::flow", "variational::flow", "flow"]}, {"name": "Rodin::Variational::Frobenius", "url": "class_rodin_1_1_variational_1_1_frobenius.html", "keys": ["rodin::variational::frobenius", "variational::frobenius", "frobenius"]}, {"name": "Rodin::Variational::FunctionBase", "url": "class_rodin_1_1_variational_1_1_function_base.html", "keys": ["rodin::variational::functionbase", "variational::functionbase", "functionbase"]}, {"name": "Rodin::Variational::FunctionBase::Cache", "url": "class_rodin_1_1_variational_1_1_function_base_1_1_cache.html", "keys": ["rodin::variational::functionbase::cache", "variational::functionbase::cache", "functionbase::cache", "cache"]}, {"name": "Rodin::Variational::FunctionBaseCopy", "url": "class_rodin_1_1_variational_1_1_function_base_copy.html", "keys": ["rodin::variational::functionbasecopy", "variational::functionbasecopy", "functionbasecopy"]}, {"name": "Rodin::Variational::GEQ", "url": "class_rodin_1_1_variational_1_1_g_e_q.html", "keys": ["rodin::variational::geq", "variational::geq", "geq"]}, {"name": "Rodin::Variational::GLL", "url": "class_rodin_1_1_variational_1_1_g_l_l.html", "keys": ["rodin::variational::gll", "variational::gll", "gll"]}, {"name": "Rodin::Variational::GLL01", "url": "class_rodin_1_1_variational_1_1_g_l_l01.html", "keys": ["rodin::variational::gll01", "variational::gll01", "gll01"]}, {"name": "Rodin::Variational::GT", "url": "class_rodin_1_1_variational_1_1_g_t.html", "keys": ["rodin::variational::gt", "variational::gt", "gt"]}, {"name": "Rodin::Variational::GlobalBilinearFormIntegratorBase", "url": "class_rodin_1_1_variational_1_1_global_bilinear_form_integrator_base.html", "keys": ["rodin::variational::globalbilinearformintegratorbase", "variational::globalbilinearformintegratorbase", "globalbilinearformintegratorbase"]}, {"name": "Rodin::Variational::Grad", "url": "class_rodin_1_1_variational_1_1_grad.html", "keys": ["rodin::variational::grad", "variational::grad", "grad"]}, {"name": "Rodin::Variational::GradBase", "url": "class_rodin_1_1_variational_1_1_grad_base.html", "keys": ["rodin::variational::gradbase", "variational::gradbase", "gradbase"]}, {"name": "Rodin::Variational::GridFunction", "url": "class_rodin_1_1_variational_1_1_grid_function.html", "keys": ["rodin::variational::gridfunction", "variational::gridfunction", "gridfunction"]}, {"name": "Rodin::Variational::GridFunctionBase", "url": "class_rodin_1_1_variational_1_1_grid_function_base.html", "keys": ["rodin::variational::gridfunctionbase", "variational::gridfunctionbase", "gridfunctionbase"]}, {"name": "Rodin::Variational::GridFunctionBaseReference", "url": "class_rodin_1_1_variational_1_1_grid_function_base_reference.html", "keys": ["rodin::variational::gridfunctionbasereference", "variational::gridfunctionbasereference", "gridfunctionbasereference"]}, {"name": "Rodin::Variational::H1", "url": "class_rodin_1_1_variational_1_1_h1.html", "keys": ["rodin::variational::h1", "variational::h1", "h1"]}, {"name": "Rodin::Variational::H1Element", "url": "class_rodin_1_1_variational_1_1_h1_element.html", "keys": ["rodin::variational::h1element", "variational::h1element", "h1element"]}, {"name": "Rodin::Variational::H1Element::BasisFunction", "url": "class_rodin_1_1_variational_1_1_h1_element_1_1_basis_function.html", "keys": ["rodin::variational::h1element::basisfunction", "variational::h1element::basisfunction", "h1element::basisfunction", "basisfunction"]}, {"name": "Rodin::Variational::H1Element::BasisFunction::DerivativeFunction", "url": "class_rodin_1_1_variational_1_1_h1_element_1_1_basis_function_1_1_derivative_function.html", "keys": ["rodin::variational::h1element::basisfunction::derivativefunction", "variational::h1element::basisfunction::derivativefunction", "h1element::basisfunction::derivativefunction", "basisfunction::derivativefunction", "derivativefunction"]}, {"name": "Rodin::Variational::H1Element::BasisFunction::GradientFunction", "url": "class_rodin_1_1_variational_1_1_h1_element_1_1_basis_function_1_1_gradient_function.html", "keys": ["rodin::variational::h1element::basisfunction::gradientfunction", "variational::h1element::basisfunction::gradientfunction", "h1element::basisfunction::gradientfunction", "basisfunction::gradientfunction", "gradientfunction"]}, {"name": "Rodin::Variational::H1Element::Evaluation", "url": "class_rodin_1_1_variational_1_1_h1_element_1_1_evaluation.html", "keys": ["rodin::variational::h1element::evaluation", "variational::h1element::evaluation", "h1element::evaluation", "evaluation"]}, {"name": "Rodin::Variational::H1Element::LinearForm", "url": "class_rodin_1_1_variational_1_1_h1_element_1_1_linear_form.html", "keys": ["rodin::variational::h1element::linearform", "variational::h1element::linearform", "h1element::linearform", "linearform"]}, {"name": "Rodin::Variational::HelmholtzForm", "url": "class_rodin_1_1_variational_1_1_helmholtz_form.html", "keys": ["rodin::variational::helmholtzform", "variational::helmholtzform", "helmholtzform"]}, {"name": "Rodin::Variational::HelmholtzForm::Kernel", "url": "class_rodin_1_1_variational_1_1_helmholtz_form_1_1_kernel.html", "keys": ["rodin::variational::helmholtzform::kernel", "variational::helmholtzform::kernel", "helmholtzform::kernel", "kernel"]}, {"name": "Rodin::Variational::IdentityMatrix", "url": "class_rodin_1_1_variational_1_1_identity_matrix.html", "keys": ["rodin::variational::identitymatrix", "variational::identitymatrix", "identitymatrix"]}, {"name": "Rodin::Variational::Im", "url": "class_rodin_1_1_variational_1_1_im.html", "keys": ["rodin::variational::im", "variational::im", "im"]}, {"name": "Rodin::Variational::IntegerFunction", "url": "class_rodin_1_1_variational_1_1_integer_function.html", "keys": ["rodin::variational::integerfunction", "variational::integerfunction", "integerfunction"]}, {"name": "Rodin::Variational::IntegerFunctionBase", "url": "class_rodin_1_1_variational_1_1_integer_function_base.html", "keys": ["rodin::variational::integerfunctionbase", "variational::integerfunctionbase", "integerfunctionbase"]}, {"name": "Rodin::Variational::Integral", "url": "class_rodin_1_1_variational_1_1_integral.html", "keys": ["rodin::variational::integral", "variational::integral", "integral"]}, {"name": "Rodin::Variational::IntegrationPoint", "url": "class_rodin_1_1_variational_1_1_integration_point.html", "keys": ["rodin::variational::integrationpoint", "variational::integrationpoint", "integrationpoint"]}, {"name": "Rodin::Variational::Integrator", "url": "class_rodin_1_1_variational_1_1_integrator.html", "keys": ["rodin::variational::integrator", "variational::integrator", "integrator"]}, {"name": "Rodin::Variational::InterfaceIntegral", "url": "class_rodin_1_1_variational_1_1_interface_integral.html", "keys": ["rodin::variational::interfaceintegral", "variational::interfaceintegral", "interfaceintegral"]}, {"name": "Rodin::Variational::InterfaceShapeFunctionBase", "url": "class_rodin_1_1_variational_1_1_interface_shape_function_base.html", "keys": ["rodin::variational::interfaceshapefunctionbase", "variational::interfaceshapefunctionbase", "interfaceshapefunctionbase"]}, {"name": "Rodin::Variational::Interpolation", "url": "class_rodin_1_1_variational_1_1_interpolation.html", "keys": ["rodin::variational::interpolation", "variational::interpolation", "interpolation"]}, {"name": "Rodin::Variational::JacobiPolynomial", "url": "class_rodin_1_1_variational_1_1_jacobi_polynomial.html", "keys": ["rodin::variational::jacobipolynomial", "variational::jacobipolynomial", "jacobipolynomial"]}, {"name": "Rodin::Variational::Jacobian", "url": "class_rodin_1_1_variational_1_1_jacobian.html", "keys": ["rodin::variational::jacobian", "variational::jacobian", "jacobian"]}, {"name": "Rodin::Variational::JacobianBase", "url": "class_rodin_1_1_variational_1_1_jacobian_base.html", "keys": ["rodin::variational::jacobianbase", "variational::jacobianbase", "jacobianbase"]}, {"name": "Rodin::Variational::Jump", "url": "class_rodin_1_1_variational_1_1_jump.html", "keys": ["rodin::variational::jump", "variational::jump", "jump"]}, {"name": "Rodin::Variational::L2", "url": "class_rodin_1_1_variational_1_1_l2.html", "keys": ["rodin::variational::l2", "variational::l2", "l2"]}, {"name": "Rodin::Variational::LEQ", "url": "class_rodin_1_1_variational_1_1_l_e_q.html", "keys": ["rodin::variational::leq", "variational::leq", "leq"]}, {"name": "Rodin::Variational::LT", "url": "class_rodin_1_1_variational_1_1_l_t.html", "keys": ["rodin::variational::lt", "variational::lt", "lt"]}, {"name": "Rodin::Variational::LagrangeBasis1D", "url": "class_rodin_1_1_variational_1_1_lagrange_basis1_d.html", "keys": ["rodin::variational::lagrangebasis1d", "variational::lagrangebasis1d", "lagrangebasis1d"]}, {"name": "Rodin::Variational::LagrangeBasisPoint", "url": "class_rodin_1_1_variational_1_1_lagrange_basis_point.html", "keys": ["rodin::variational::lagrangebasispoint", "variational::lagrangebasispoint", "lagrangebasispoint"]}, {"name": "Rodin::Variational::LagrangeBasisQuadrilateral", "url": "class_rodin_1_1_variational_1_1_lagrange_basis_quadrilateral.html", "keys": ["rodin::variational::lagrangebasisquadrilateral", "variational::lagrangebasisquadrilateral", "lagrangebasisquadrilateral"]}, {"name": "Rodin::Variational::LagrangeBasisSegment", "url": "class_rodin_1_1_variational_1_1_lagrange_basis_segment.html", "keys": ["rodin::variational::lagrangebasissegment", "variational::lagrangebasissegment", "lagrangebasissegment"]}, {"name": "Rodin::Variational::LagrangeBasisTetrahedron", "url": "class_rodin_1_1_variational_1_1_lagrange_basis_tetrahedron.html", "keys": ["rodin::variational::lagrangebasistetrahedron", "variational::lagrangebasistetrahedron", "lagrangebasistetrahedron"]}, {"name": "Rodin::Variational::LagrangeBasisTriangle", "url": "class_rodin_1_1_variational_1_1_lagrange_basis_triangle.html", "keys": ["rodin::variational::lagrangebasistriangle", "variational::lagrangebasistriangle", "lagrangebasistriangle"]}, {"name": "Rodin::Variational::LagrangeBasisWedge", "url": "class_rodin_1_1_variational_1_1_lagrange_basis_wedge.html", "keys": ["rodin::variational::lagrangebasiswedge", "variational::lagrangebasiswedge", "lagrangebasiswedge"]}, {"name": "Rodin::Variational::LagrangeFiniteElementBase", "url": "class_rodin_1_1_variational_1_1_lagrange_finite_element_base.html", "keys": ["rodin::variational::lagrangefiniteelementbase", "variational::lagrangefiniteelementbase", "lagrangefiniteelementbase"]}, {"name": "Rodin::Variational::LegendrePolynomial", "url": "class_rodin_1_1_variational_1_1_legendre_polynomial.html", "keys": ["rodin::variational::legendrepolynomial", "variational::legendrepolynomial", "legendrepolynomial"]}, {"name": "Rodin::Variational::LinearElasticityForm", "url": "class_rodin_1_1_variational_1_1_linear_elasticity_form.html", "keys": ["rodin::variational::linearelasticityform", "variational::linearelasticityform", "linearelasticityform"]}, {"name": "Rodin::Variational::LinearElasticityForm::Kernel", "url": "class_rodin_1_1_variational_1_1_linear_elasticity_form_1_1_kernel.html", "keys": ["rodin::variational::linearelasticityform::kernel", "variational::linearelasticityform::kernel", "linearelasticityform::kernel", "kernel"]}, {"name": "Rodin::Variational::LinearElasticityIntegral", "url": "class_rodin_1_1_variational_1_1_linear_elasticity_integral.html", "keys": ["rodin::variational::linearelasticityintegral", "variational::linearelasticityintegral", "linearelasticityintegral"]}, {"name": "Rodin::Variational::LinearElasticityIntegrator", "url": "class_rodin_1_1_variational_1_1_linear_elasticity_integrator.html", "keys": ["rodin::variational::linearelasticityintegrator", "variational::linearelasticityintegrator", "linearelasticityintegrator"]}, {"name": "Rodin::Variational::LinearForm", "url": "class_rodin_1_1_variational_1_1_linear_form.html", "keys": ["rodin::variational::linearform", "variational::linearform", "linearform"]}, {"name": "Rodin::Variational::LinearFormBase", "url": "class_rodin_1_1_variational_1_1_linear_form_base.html", "keys": ["rodin::variational::linearformbase", "variational::linearformbase", "linearformbase"]}, {"name": "Rodin::Variational::LinearFormIntegratorBase", "url": "class_rodin_1_1_variational_1_1_linear_form_integrator_base.html", "keys": ["rodin::variational::linearformintegratorbase", "variational::linearformintegratorbase", "linearformintegratorbase"]}, {"name": "Rodin::Variational::LocalBilinearFormIntegratorBase", "url": "class_rodin_1_1_variational_1_1_local_bilinear_form_integrator_base.html", "keys": ["rodin::variational::localbilinearformintegratorbase", "variational::localbilinearformintegratorbase", "localbilinearformintegratorbase"]}, {"name": "Rodin::Variational::MassForm", "url": "class_rodin_1_1_variational_1_1_mass_form.html", "keys": ["rodin::variational::massform", "variational::massform", "massform"]}, {"name": "Rodin::Variational::MassForm::Kernel", "url": "class_rodin_1_1_variational_1_1_mass_form_1_1_kernel.html", "keys": ["rodin::variational::massform::kernel", "variational::massform::kernel", "massform::kernel", "kernel"]}, {"name": "Rodin::Variational::MatrixFunction", "url": "class_rodin_1_1_variational_1_1_matrix_function.html", "keys": ["rodin::variational::matrixfunction", "variational::matrixfunction", "matrixfunction"]}, {"name": "Rodin::Variational::MatrixFunctionBase", "url": "class_rodin_1_1_variational_1_1_matrix_function_base.html", "keys": ["rodin::variational::matrixfunctionbase", "variational::matrixfunctionbase", "matrixfunctionbase"]}, {"name": "Rodin::Variational::Max", "url": "class_rodin_1_1_variational_1_1_max.html", "keys": ["rodin::variational::max", "variational::max", "max"]}, {"name": "Rodin::Variational::Min", "url": "class_rodin_1_1_variational_1_1_min.html", "keys": ["rodin::variational::min", "variational::min", "min"]}, {"name": "Rodin::Variational::Mult", "url": "class_rodin_1_1_variational_1_1_mult.html", "keys": ["rodin::variational::mult", "variational::mult", "mult"]}, {"name": "Rodin::Variational::NEQ", "url": "class_rodin_1_1_variational_1_1_n_e_q.html", "keys": ["rodin::variational::neq", "variational::neq", "neq"]}, {"name": "Rodin::Variational::NamedFormStorage", "url": "class_rodin_1_1_variational_1_1_named_form_storage.html", "keys": ["rodin::variational::namedformstorage", "variational::namedformstorage", "namedformstorage"]}, {"name": "Rodin::Variational::NonLinearFormIntegrator", "url": "class_rodin_1_1_variational_1_1_non_linear_form_integrator.html", "keys": ["rodin::variational::nonlinearformintegrator", "variational::nonlinearformintegrator", "nonlinearformintegrator"]}, {"name": "Rodin::Variational::OR", "url": "class_rodin_1_1_variational_1_1_o_r.html", "keys": ["rodin::variational::or", "variational::or", "or"]}, {"name": "Rodin::Variational::P0", "url": "class_rodin_1_1_variational_1_1_p0.html", "keys": ["rodin::variational::p0", "variational::p0", "p0"]}, {"name": "Rodin::Variational::P0Element", "url": "class_rodin_1_1_variational_1_1_p0_element.html", "keys": ["rodin::variational::p0element", "variational::p0element", "p0element"]}, {"name": "Rodin::Variational::P0Element::BasisFunction", "url": "class_rodin_1_1_variational_1_1_p0_element_1_1_basis_function.html", "keys": ["rodin::variational::p0element::basisfunction", "variational::p0element::basisfunction", "p0element::basisfunction", "basisfunction"]}, {"name": "Rodin::Variational::P0Element::BasisFunction::DerivativeFunction", "url": "class_rodin_1_1_variational_1_1_p0_element_1_1_basis_function_1_1_derivative_function.html", "keys": ["rodin::variational::p0element::basisfunction::derivativefunction", "variational::p0element::basisfunction::derivativefunction", "p0element::basisfunction::derivativefunction", "basisfunction::derivativefunction", "derivativefunction"]}, {"name": "Rodin::Variational::P0Element::LinearForm", "url": "class_rodin_1_1_variational_1_1_p0_element_1_1_linear_form.html", "keys": ["rodin::variational::p0element::linearform", "variational::p0element::linearform", "p0element::linearform", "linearform"]}, {"name": "Rodin::Variational::P0g", "url": "class_rodin_1_1_variational_1_1_p0g.html", "keys": ["rodin::variational::p0g", "variational::p0g", "p0g"]}, {"name": "Rodin::Variational::P0gElement", "url": "class_rodin_1_1_variational_1_1_p0g_element.html", "keys": ["rodin::variational::p0gelement", "variational::p0gelement", "p0gelement"]}, {"name": "Rodin::Variational::P0gElement::BasisFunction", "url": "class_rodin_1_1_variational_1_1_p0g_element_1_1_basis_function.html", "keys": ["rodin::variational::p0gelement::basisfunction", "variational::p0gelement::basisfunction", "p0gelement::basisfunction", "basisfunction"]}, {"name": "Rodin::Variational::P0gElement::BasisFunction::DerivativeFunction", "url": "class_rodin_1_1_variational_1_1_p0g_element_1_1_basis_function_1_1_derivative_function.html", "keys": ["rodin::variational::p0gelement::basisfunction::derivativefunction", "variational::p0gelement::basisfunction::derivativefunction", "p0gelement::basisfunction::derivativefunction", "basisfunction::derivativefunction", "derivativefunction"]}, {"name": "Rodin::Variational::P0gElement::LinearForm", "url": "class_rodin_1_1_variational_1_1_p0g_element_1_1_linear_form.html", "keys": ["rodin::variational::p0gelement::linearform", "variational::p0gelement::linearform", "p0gelement::linearform", "linearform"]}, {"name": "Rodin::Variational::P1", "url": "class_rodin_1_1_variational_1_1_p1.html", "keys": ["rodin::variational::p1", "variational::p1", "p1"]}, {"name": "Rodin::Variational::P1Element", "url": "class_rodin_1_1_variational_1_1_p1_element.html", "keys": ["rodin::variational::p1element", "variational::p1element", "p1element"]}, {"name": "Rodin::Variational::P1Element::BasisFunction", "url": "class_rodin_1_1_variational_1_1_p1_element_1_1_basis_function.html", "keys": ["rodin::variational::p1element::basisfunction", "variational::p1element::basisfunction", "p1element::basisfunction", "basisfunction"]}, {"name": "Rodin::Variational::P1Element::BasisFunction::DerivativeFunction", "url": "class_rodin_1_1_variational_1_1_p1_element_1_1_basis_function_1_1_derivative_function.html", "keys": ["rodin::variational::p1element::basisfunction::derivativefunction", "variational::p1element::basisfunction::derivativefunction", "p1element::basisfunction::derivativefunction", "basisfunction::derivativefunction", "derivativefunction"]}, {"name": "Rodin::Variational::P1Element::BasisFunction::GradientFunction", "url": "class_rodin_1_1_variational_1_1_p1_element_1_1_basis_function_1_1_gradient_function.html", "keys": ["rodin::variational::p1element::basisfunction::gradientfunction", "variational::p1element::basisfunction::gradientfunction", "p1element::basisfunction::gradientfunction", "basisfunction::gradientfunction", "gradientfunction"]}, {"name": "Rodin::Variational::P1Element::LinearForm", "url": "class_rodin_1_1_variational_1_1_p1_element_1_1_linear_form.html", "keys": ["rodin::variational::p1element::linearform", "variational::p1element::linearform", "p1element::linearform", "linearform"]}, {"name": "Rodin::Variational::PeriodicBC", "url": "class_rodin_1_1_variational_1_1_periodic_b_c.html", "keys": ["rodin::variational::periodicbc", "variational::periodicbc", "periodicbc"]}, {"name": "Rodin::Variational::PeriodicBCBase", "url": "class_rodin_1_1_variational_1_1_periodic_b_c_base.html", "keys": ["rodin::variational::periodicbcbase", "variational::periodicbcbase", "periodicbcbase"]}, {"name": "Rodin::Variational::Potential", "url": "class_rodin_1_1_variational_1_1_potential.html", "keys": ["rodin::variational::potential", "variational::potential", "potential"]}, {"name": "Rodin::Variational::Pow", "url": "class_rodin_1_1_variational_1_1_pow.html", "keys": ["rodin::variational::pow", "variational::pow", "pow"]}, {"name": "Rodin::Variational::Problem", "url": "class_rodin_1_1_variational_1_1_problem.html", "keys": ["rodin::variational::problem", "variational::problem", "problem"]}, {"name": "Rodin::Variational::ProblemBase", "url": "class_rodin_1_1_variational_1_1_problem_base.html", "keys": ["rodin::variational::problembase", "variational::problembase", "problembase"]}, {"name": "Rodin::Variational::ProblemBody", "url": "class_rodin_1_1_variational_1_1_problem_body.html", "keys": ["rodin::variational::problembody", "variational::problembody", "problembody"]}, {"name": "Rodin::Variational::ProblemBodyBase", "url": "class_rodin_1_1_variational_1_1_problem_body_base.html", "keys": ["rodin::variational::problembodybase", "variational::problembodybase", "problembodybase"]}, {"name": "Rodin::Variational::ProblemUVBase", "url": "class_rodin_1_1_variational_1_1_problem_u_v_base.html", "keys": ["rodin::variational::problemuvbase", "variational::problemuvbase", "problemuvbase"]}, {"name": "Rodin::Variational::ProblemUsBase", "url": "class_rodin_1_1_variational_1_1_problem_us_base.html", "keys": ["rodin::variational::problemusbase", "variational::problemusbase", "problemusbase"]}, {"name": "Rodin::Variational::PyramidModal", "url": "class_rodin_1_1_variational_1_1_pyramid_modal.html", "keys": ["rodin::variational::pyramidmodal", "variational::pyramidmodal", "pyramidmodal"]}, {"name": "Rodin::Variational::QuadratureRule", "url": "class_rodin_1_1_variational_1_1_quadrature_rule.html", "keys": ["rodin::variational::quadraturerule", "variational::quadraturerule", "quadraturerule"]}, {"name": "Rodin::Variational::Re", "url": "class_rodin_1_1_variational_1_1_re.html", "keys": ["rodin::variational::re", "variational::re", "re"]}, {"name": "Rodin::Variational::RealFunction", "url": "class_rodin_1_1_variational_1_1_real_function.html", "keys": ["rodin::variational::realfunction", "variational::realfunction", "realfunction"]}, {"name": "Rodin::Variational::RealFunctionBase", "url": "class_rodin_1_1_variational_1_1_real_function_base.html", "keys": ["rodin::variational::realfunctionbase", "variational::realfunctionbase", "realfunctionbase"]}, {"name": "Rodin::Variational::ReferenceTabulation", "url": "class_rodin_1_1_variational_1_1_reference_tabulation.html", "keys": ["rodin::variational::referencetabulation", "variational::referencetabulation", "referencetabulation"]}, {"name": "Rodin::Variational::ReferenceTabulationCache", "url": "class_rodin_1_1_variational_1_1_reference_tabulation_cache.html", "keys": ["rodin::variational::referencetabulationcache", "variational::referencetabulationcache", "referencetabulationcache"]}, {"name": "Rodin::Variational::RelativeError", "url": "class_rodin_1_1_variational_1_1_relative_error.html", "keys": ["rodin::variational::relativeerror", "variational::relativeerror", "relativeerror"]}, {"name": "Rodin::Variational::ScalarFunctionBase", "url": "class_rodin_1_1_variational_1_1_scalar_function_base.html", "keys": ["rodin::variational::scalarfunctionbase", "variational::scalarfunctionbase", "scalarfunctionbase"]}, {"name": "Rodin::Variational::ShapeFunction", "url": "class_rodin_1_1_variational_1_1_shape_function.html", "keys": ["rodin::variational::shapefunction", "variational::shapefunction", "shapefunction"]}, {"name": "Rodin::Variational::ShapeFunctionBase", "url": "class_rodin_1_1_variational_1_1_shape_function_base.html", "keys": ["rodin::variational::shapefunctionbase", "variational::shapefunctionbase", "shapefunctionbase"]}, {"name": "Rodin::Variational::Sin", "url": "class_rodin_1_1_variational_1_1_sin.html", "keys": ["rodin::variational::sin", "variational::sin", "sin"]}, {"name": "Rodin::Variational::Sinh", "url": "class_rodin_1_1_variational_1_1_sinh.html", "keys": ["rodin::variational::sinh", "variational::sinh", "sinh"]}, {"name": "Rodin::Variational::SparseProblem", "url": "class_rodin_1_1_variational_1_1_sparse_problem.html", "keys": ["rodin::variational::sparseproblem", "variational::sparseproblem", "sparseproblem"]}, {"name": "Rodin::Variational::Sqrt", "url": "class_rodin_1_1_variational_1_1_sqrt.html", "keys": ["rodin::variational::sqrt", "variational::sqrt", "sqrt"]}, {"name": "Rodin::Variational::Sum", "url": "class_rodin_1_1_variational_1_1_sum.html", "keys": ["rodin::variational::sum", "variational::sum", "sum"]}, {"name": "Rodin::Variational::Tan", "url": "class_rodin_1_1_variational_1_1_tan.html", "keys": ["rodin::variational::tan", "variational::tan", "tan"]}, {"name": "Rodin::Variational::TensorFunction", "url": "class_rodin_1_1_variational_1_1_tensor_function.html", "keys": ["rodin::variational::tensorfunction", "variational::tensorfunction", "tensorfunction"]}, {"name": "Rodin::Variational::TestFunction", "url": "class_rodin_1_1_variational_1_1_test_function.html", "keys": ["rodin::variational::testfunction", "variational::testfunction", "testfunction"]}, {"name": "Rodin::Variational::TestFunctionMismatchException", "url": "class_rodin_1_1_variational_1_1_test_function_mismatch_exception.html", "keys": ["rodin::variational::testfunctionmismatchexception", "variational::testfunctionmismatchexception", "testfunctionmismatchexception"]}, {"name": "Rodin::Variational::TetrahedronBlend", "url": "class_rodin_1_1_variational_1_1_tetrahedron_blend.html", "keys": ["rodin::variational::tetrahedronblend", "variational::tetrahedronblend", "tetrahedronblend"]}, {"name": "Rodin::Variational::Trace", "url": "class_rodin_1_1_variational_1_1_trace.html", "keys": ["rodin::variational::trace", "variational::trace", "trace"]}, {"name": "Rodin::Variational::TraceOperator", "url": "class_rodin_1_1_variational_1_1_trace_operator.html", "keys": ["rodin::variational::traceoperator", "variational::traceoperator", "traceoperator"]}, {"name": "Rodin::Variational::Transpose", "url": "class_rodin_1_1_variational_1_1_transpose.html", "keys": ["rodin::variational::transpose", "variational::transpose", "transpose"]}, {"name": "Rodin::Variational::TrialFunction", "url": "class_rodin_1_1_variational_1_1_trial_function.html", "keys": ["rodin::variational::trialfunction", "variational::trialfunction", "trialfunction"]}, {"name": "Rodin::Variational::TrialFunctionMismatchException", "url": "class_rodin_1_1_variational_1_1_trial_function_mismatch_exception.html", "keys": ["rodin::variational::trialfunctionmismatchexception", "variational::trialfunctionmismatchexception", "trialfunctionmismatchexception"]}, {"name": "Rodin::Variational::TrialFunctionReference", "url": "class_rodin_1_1_variational_1_1_trial_function_reference.html", "keys": ["rodin::variational::trialfunctionreference", "variational::trialfunctionreference", "trialfunctionreference"]}, {"name": "Rodin::Variational::TriangleBlend", "url": "class_rodin_1_1_variational_1_1_triangle_blend.html", "keys": ["rodin::variational::triangleblend", "variational::triangleblend", "triangleblend"]}, {"name": "Rodin::Variational::UnaryMinus", "url": "class_rodin_1_1_variational_1_1_unary_minus.html", "keys": ["rodin::variational::unaryminus", "variational::unaryminus", "unaryminus"]}, {"name": "Rodin::Variational::UndeterminedTraceDomainException", "url": "class_rodin_1_1_variational_1_1_undetermined_trace_domain_exception.html", "keys": ["rodin::variational::undeterminedtracedomainexception", "variational::undeterminedtracedomainexception", "undeterminedtracedomainexception"]}, {"name": "Rodin::Variational::VandermondePyramid", "url": "class_rodin_1_1_variational_1_1_vandermonde_pyramid.html", "keys": ["rodin::variational::vandermondepyramid", "variational::vandermondepyramid", "vandermondepyramid"]}, {"name": "Rodin::Variational::VandermondeTetrahedron", "url": "class_rodin_1_1_variational_1_1_vandermonde_tetrahedron.html", "keys": ["rodin::variational::vandermondetetrahedron", "variational::vandermondetetrahedron", "vandermondetetrahedron"]}, {"name": "Rodin::Variational::VandermondeTriangle", "url": "class_rodin_1_1_variational_1_1_vandermonde_triangle.html", "keys": ["rodin::variational::vandermondetriangle", "variational::vandermondetriangle", "vandermondetriangle"]}, {"name": "Rodin::Variational::VectorFunction", "url": "class_rodin_1_1_variational_1_1_vector_function.html", "keys": ["rodin::variational::vectorfunction", "variational::vectorfunction", "vectorfunction"]}, {"name": "Rodin::Variational::VectorFunctionBase", "url": "class_rodin_1_1_variational_1_1_vector_function_base.html", "keys": ["rodin::variational::vectorfunctionbase", "variational::vectorfunctionbase", "vectorfunctionbase"]}, {"name": "Rodin::Variational::WarpBlendTetrahedron", "url": "class_rodin_1_1_variational_1_1_warp_blend_tetrahedron.html", "keys": ["rodin::variational::warpblendtetrahedron", "variational::warpblendtetrahedron", "warpblendtetrahedron"]}, {"name": "Rodin::Variational::WarpBlendTriangle", "url": "class_rodin_1_1_variational_1_1_warp_blend_triangle.html", "keys": ["rodin::variational::warpblendtriangle", "variational::warpblendtriangle", "warpblendtriangle"]}, {"name": "Rodin::Variational::WarpFactor1D", "url": "class_rodin_1_1_variational_1_1_warp_factor1_d.html", "keys": ["rodin::variational::warpfactor1d", "variational::warpfactor1d", "warpfactor1d"]}, {"name": "Rodin::Variational::WarpShiftFace2D", "url": "class_rodin_1_1_variational_1_1_warp_shift_face2_d.html", "keys": ["rodin::variational::warpshiftface2d", "variational::warpshiftface2d", "warpshiftface2d"]}, {"name": "Rodin::Variational::WarpShiftFace3D", "url": "class_rodin_1_1_variational_1_1_warp_shift_face3_d.html", "keys": ["rodin::variational::warpshiftface3d", "variational::warpshiftface3d", "warpshiftface3d"]}, {"name": "Rodin::Variational::Zero", "url": "class_rodin_1_1_variational_1_1_zero.html", "keys": ["rodin::variational::zero", "variational::zero", "zero"]}, {"name": "Rodin::Adaptation::SWIFT::AdmissibilityReport", "url": "struct_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_admissibility_report.html", "keys": ["rodin::adaptation::swift::admissibilityreport", "adaptation::swift::admissibilityreport", "swift::admissibilityreport", "admissibilityreport"]}, {"name": "Rodin::Adaptation::SWIFT::Parameters", "url": "struct_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_parameters.html", "keys": ["rodin::adaptation::swift::parameters", "adaptation::swift::parameters", "swift::parameters", "parameters"]}, {"name": "Rodin::Adaptation::SWIFT::Parameters::Convergence", "url": "struct_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_parameters_1_1_convergence.html", "keys": ["rodin::adaptation::swift::parameters::convergence", "adaptation::swift::parameters::convergence", "swift::parameters::convergence", "parameters::convergence", "convergence"]}, {"name": "Rodin::Adaptation::SWIFT::Parameters::Convergence::Iterations", "url": "struct_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_parameters_1_1_convergence_1_1_iterations.html", "keys": ["rodin::adaptation::swift::parameters::convergence::iterations", "adaptation::swift::parameters::convergence::iterations", "swift::parameters::convergence::iterations", "parameters::convergence::iterations", "convergence::iterations", "iterations"]}, {"name": "Rodin::Adaptation::SWIFT::Parameters::Convergence::Tolerance", "url": "struct_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_parameters_1_1_convergence_1_1_tolerance.html", "keys": ["rodin::adaptation::swift::parameters::convergence::tolerance", "adaptation::swift::parameters::convergence::tolerance", "swift::parameters::convergence::tolerance", "parameters::convergence::tolerance", "convergence::tolerance", "tolerance"]}, {"name": "Rodin::Adaptation::SWIFT::Parameters::Globalization", "url": "struct_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_parameters_1_1_globalization.html", "keys": ["rodin::adaptation::swift::parameters::globalization", "adaptation::swift::parameters::globalization", "swift::parameters::globalization", "parameters::globalization", "globalization"]}, {"name": "Rodin::Adaptation::SWIFT::Parameters::Linear", "url": "struct_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_parameters_1_1_linear.html", "keys": ["rodin::adaptation::swift::parameters::linear", "adaptation::swift::parameters::linear", "swift::parameters::linear", "parameters::linear", "linear"]}, {"name": "Rodin::Adaptation::SWIFT::Parameters::Model", "url": "struct_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_parameters_1_1_model.html", "keys": ["rodin::adaptation::swift::parameters::model", "adaptation::swift::parameters::model", "swift::parameters::model", "parameters::model", "model"]}, {"name": "Rodin::Adaptation::SWIFT::Parameters::Model::Distribution", "url": "struct_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_parameters_1_1_model_1_1_distribution.html", "keys": ["rodin::adaptation::swift::parameters::model::distribution", "adaptation::swift::parameters::model::distribution", "swift::parameters::model::distribution", "parameters::model::distribution", "model::distribution", "distribution"]}, {"name": "Rodin::Adaptation::SWIFT::Parameters::Model::Hinge", "url": "struct_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_parameters_1_1_model_1_1_hinge.html", "keys": ["rodin::adaptation::swift::parameters::model::hinge", "adaptation::swift::parameters::model::hinge", "swift::parameters::model::hinge", "parameters::model::hinge", "model::hinge", "hinge"]}, {"name": "Rodin::Adaptation::SWIFT::Parameters::Quadrature", "url": "struct_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_parameters_1_1_quadrature.html", "keys": ["rodin::adaptation::swift::parameters::quadrature", "adaptation::swift::parameters::quadrature", "swift::parameters::quadrature", "parameters::quadrature", "quadrature"]}, {"name": "Rodin::Adaptation::SWIFT::Parameters::Sampling", "url": "struct_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_parameters_1_1_sampling.html", "keys": ["rodin::adaptation::swift::parameters::sampling", "adaptation::swift::parameters::sampling", "swift::parameters::sampling", "parameters::sampling", "sampling"]}, {"name": "Rodin::Adaptation::SWIFT::Report", "url": "struct_rodin_1_1_adaptation_1_1_s_w_i_f_t_1_1_report.html", "keys": ["rodin::adaptation::swift::report", "adaptation::swift::report", "swift::report", "report"]}, {"name": "Rodin::Alert::BlueT", "url": "struct_rodin_1_1_alert_1_1_blue_t.html", "keys": ["rodin::alert::bluet", "alert::bluet", "bluet"]}, {"name": "Rodin::Alert::BoldT", "url": "struct_rodin_1_1_alert_1_1_bold_t.html", "keys": ["rodin::alert::boldt", "alert::boldt", "boldt"]}, {"name": "Rodin::Alert::BrightGrayT", "url": "struct_rodin_1_1_alert_1_1_bright_gray_t.html", "keys": ["rodin::alert::brightgrayt", "alert::brightgrayt", "brightgrayt"]}, {"name": "Rodin::Alert::BrightWhiteT", "url": "struct_rodin_1_1_alert_1_1_bright_white_t.html", "keys": ["rodin::alert::brightwhitet", "alert::brightwhitet", "brightwhitet"]}, {"name": "Rodin::Alert::CyanT", "url": "struct_rodin_1_1_alert_1_1_cyan_t.html", "keys": ["rodin::alert::cyant", "alert::cyant", "cyant"]}, {"name": "Rodin::Alert::GrayT", "url": "struct_rodin_1_1_alert_1_1_gray_t.html", "keys": ["rodin::alert::grayt", "alert::grayt", "grayt"]}, {"name": "Rodin::Alert::GreenT", "url": "struct_rodin_1_1_alert_1_1_green_t.html", "keys": ["rodin::alert::greent", "alert::greent", "greent"]}, {"name": "Rodin::Alert::MagentaT", "url": "struct_rodin_1_1_alert_1_1_magenta_t.html", "keys": ["rodin::alert::magentat", "alert::magentat", "magentat"]}, {"name": "Rodin::Alert::NewLineT", "url": "struct_rodin_1_1_alert_1_1_new_line_t.html", "keys": ["rodin::alert::newlinet", "alert::newlinet", "newlinet"]}, {"name": "Rodin::Alert::NoColorT", "url": "struct_rodin_1_1_alert_1_1_no_color_t.html", "keys": ["rodin::alert::nocolort", "alert::nocolort", "nocolort"]}, {"name": "Rodin::Alert::RGB", "url": "struct_rodin_1_1_alert_1_1_r_g_b.html", "keys": ["rodin::alert::rgb", "alert::rgb", "rgb"]}, {"name": "Rodin::Alert::RaiseT", "url": "struct_rodin_1_1_alert_1_1_raise_t.html", "keys": ["rodin::alert::raiset", "alert::raiset", "raiset"]}, {"name": "Rodin::Alert::RedT", "url": "struct_rodin_1_1_alert_1_1_red_t.html", "keys": ["rodin::alert::redt", "alert::redt", "redt"]}, {"name": "Rodin::Alert::ResetT", "url": "struct_rodin_1_1_alert_1_1_reset_t.html", "keys": ["rodin::alert::resett", "alert::resett", "resett"]}, {"name": "Rodin::Alert::StylizeT", "url": "struct_rodin_1_1_alert_1_1_stylize_t.html", "keys": ["rodin::alert::stylizet", "alert::stylizet", "stylizet"]}, {"name": "Rodin::Alert::WhiteT", "url": "struct_rodin_1_1_alert_1_1_white_t.html", "keys": ["rodin::alert::whitet", "alert::whitet", "whitet"]}, {"name": "Rodin::Alert::YellowT", "url": "struct_rodin_1_1_alert_1_1_yellow_t.html", "keys": ["rodin::alert::yellowt", "alert::yellowt", "yellowt"]}, {"name": "Rodin::Assembly::ConstraintMap::Entry", "url": "struct_rodin_1_1_assembly_1_1_constraint_map_1_1_entry.html", "keys": ["rodin::assembly::constraintmap::entry", "assembly::constraintmap::entry", "constraintmap::entry", "entry"]}, {"name": "Rodin::DefaultSentinelT", "url": "struct_rodin_1_1_default_sentinel_t.html", "keys": ["rodin::defaultsentinelt", "defaultsentinelt"]}, {"name": "Rodin::FormLanguage::ColsAtCompileTime", "url": "struct_rodin_1_1_form_language_1_1_cols_at_compile_time.html", "keys": ["rodin::formlanguage::colsatcompiletime", "formlanguage::colsatcompiletime", "colsatcompiletime"]}, {"name": "Rodin::FormLanguage::Debug", "url": "struct_rodin_1_1_form_language_1_1_debug.html", "keys": ["rodin::formlanguage::debug", "formlanguage::debug", "debug"]}, {"name": "Rodin::FormLanguage::Division", "url": "struct_rodin_1_1_form_language_1_1_division.html", "keys": ["rodin::formlanguage::division", "formlanguage::division", "division"]}, {"name": "Rodin::FormLanguage::Dot", "url": "struct_rodin_1_1_form_language_1_1_dot.html", "keys": ["rodin::formlanguage::dot", "formlanguage::dot", "dot"]}, {"name": "Rodin::FormLanguage::FunctionDerived", "url": "struct_rodin_1_1_form_language_1_1_function_derived.html", "keys": ["rodin::formlanguage::functionderived", "formlanguage::functionderived", "functionderived"]}, {"name": "Rodin::FormLanguage::IsEigenObject", "url": "struct_rodin_1_1_form_language_1_1_is_eigen_object.html", "keys": ["rodin::formlanguage::iseigenobject", "formlanguage::iseigenobject", "iseigenobject"]}, {"name": "Rodin::FormLanguage::IsInterfaceShapeFunction", "url": "struct_rodin_1_1_form_language_1_1_is_interface_shape_function.html", "keys": ["rodin::formlanguage::isinterfaceshapefunction", "formlanguage::isinterfaceshapefunction", "isinterfaceshapefunction"]}, {"name": "Rodin::FormLanguage::IsMatrixRange", "url": "struct_rodin_1_1_form_language_1_1_is_matrix_range.html", "keys": ["rodin::formlanguage::ismatrixrange", "formlanguage::ismatrixrange", "ismatrixrange"]}, {"name": "Rodin::FormLanguage::IsNamedForm", "url": "struct_rodin_1_1_form_language_1_1_is_named_form.html", "keys": ["rodin::formlanguage::isnamedform", "formlanguage::isnamedform", "isnamedform"]}, {"name": "Rodin::FormLanguage::IsPlainObject", "url": "struct_rodin_1_1_form_language_1_1_is_plain_object.html", "keys": ["rodin::formlanguage::isplainobject", "formlanguage::isplainobject", "isplainobject"]}, {"name": "Rodin::FormLanguage::IsSpatialMatrix", "url": "struct_rodin_1_1_form_language_1_1_is_spatial_matrix.html", "keys": ["rodin::formlanguage::isspatialmatrix", "formlanguage::isspatialmatrix", "isspatialmatrix"]}, {"name": "Rodin::FormLanguage::IsSpatialVector", "url": "struct_rodin_1_1_form_language_1_1_is_spatial_vector.html", "keys": ["rodin::formlanguage::isspatialvector", "formlanguage::isspatialvector", "isspatialvector"]}, {"name": "Rodin::FormLanguage::IsSpecialized", "url": "struct_rodin_1_1_form_language_1_1_is_specialized.html", "keys": ["rodin::formlanguage::isspecialized", "formlanguage::isspecialized", "isspecialized"]}, {"name": "Rodin::FormLanguage::IsTensorRange", "url": "struct_rodin_1_1_form_language_1_1_is_tensor_range.html", "keys": ["rodin::formlanguage::istensorrange", "formlanguage::istensorrange", "istensorrange"]}, {"name": "Rodin::FormLanguage::IsVectorRange", "url": "struct_rodin_1_1_form_language_1_1_is_vector_range.html", "keys": ["rodin::formlanguage::isvectorrange", "formlanguage::isvectorrange", "isvectorrange"]}, {"name": "Rodin::FormLanguage::Minus", "url": "struct_rodin_1_1_form_language_1_1_minus.html", "keys": ["rodin::formlanguage::minus", "formlanguage::minus", "minus"]}, {"name": "Rodin::FormLanguage::Mult", "url": "struct_rodin_1_1_form_language_1_1_mult.html", "keys": ["rodin::formlanguage::mult", "formlanguage::mult", "mult"]}, {"name": "Rodin::FormLanguage::NamedFormOperatorType", "url": "struct_rodin_1_1_form_language_1_1_named_form_operator_type.html", "keys": ["rodin::formlanguage::namedformoperatortype", "formlanguage::namedformoperatortype", "namedformoperatortype"]}, {"name": "Rodin::FormLanguage::RangeKindOf", "url": "struct_rodin_1_1_form_language_1_1_range_kind_of.html", "keys": ["rodin::formlanguage::rangekindof", "formlanguage::rangekindof", "rangekindof"]}, {"name": "Rodin::FormLanguage::RangeOf", "url": "struct_rodin_1_1_form_language_1_1_range_of.html", "keys": ["rodin::formlanguage::rangeof", "formlanguage::rangeof", "rangeof"]}, {"name": "Rodin::FormLanguage::ResultOf", "url": "struct_rodin_1_1_form_language_1_1_result_of.html", "keys": ["rodin::formlanguage::resultof", "formlanguage::resultof", "resultof"]}, {"name": "Rodin::FormLanguage::Sum", "url": "struct_rodin_1_1_form_language_1_1_sum.html", "keys": ["rodin::formlanguage::sum", "formlanguage::sum", "sum"]}, {"name": "Rodin::FormLanguage::Traits", "url": "struct_rodin_1_1_form_language_1_1_traits.html", "keys": ["rodin::formlanguage::traits", "formlanguage::traits", "traits"]}, {"name": "Rodin::FormLanguage::UnaryMinus", "url": "struct_rodin_1_1_form_language_1_1_unary_minus.html", "keys": ["rodin::formlanguage::unaryminus", "formlanguage::unaryminus", "unaryminus"]}, {"name": "Rodin::Geometry::Euclidean::HasConnectMethod", "url": "struct_rodin_1_1_geometry_1_1_euclidean_1_1_has_connect_method.html", "keys": ["rodin::geometry::euclidean::hasconnectmethod", "geometry::euclidean::hasconnectmethod", "euclidean::hasconnectmethod", "hasconnectmethod"]}, {"name": "Rodin::Geometry::Euclidean::HasDistanceMethod", "url": "struct_rodin_1_1_geometry_1_1_euclidean_1_1_has_distance_method.html", "keys": ["rodin::geometry::euclidean::hasdistancemethod", "geometry::euclidean::hasdistancemethod", "euclidean::hasdistancemethod", "hasdistancemethod"]}, {"name": "Rodin::Geometry::Euclidean::HasIntersectMethod", "url": "struct_rodin_1_1_geometry_1_1_euclidean_1_1_has_intersect_method.html", "keys": ["rodin::geometry::euclidean::hasintersectmethod", "geometry::euclidean::hasintersectmethod", "euclidean::hasintersectmethod", "hasintersectmethod"]}, {"name": "Rodin::Geometry::Polytope::Key::SymmetricEquality", "url": "struct_rodin_1_1_geometry_1_1_polytope_1_1_key_1_1_symmetric_equality.html", "keys": ["rodin::geometry::polytope::key::symmetricequality", "geometry::polytope::key::symmetricequality", "polytope::key::symmetricequality", "key::symmetricequality", "symmetricequality"]}, {"name": "Rodin::Geometry::Polytope::Key::SymmetricHash", "url": "struct_rodin_1_1_geometry_1_1_polytope_1_1_key_1_1_symmetric_hash.html", "keys": ["rodin::geometry::polytope::key::symmetrichash", "geometry::polytope::key::symmetrichash", "polytope::key::symmetrichash", "key::symmetrichash", "symmetrichash"]}, {"name": "Rodin::Geometry::Polytope::Project", "url": "struct_rodin_1_1_geometry_1_1_polytope_1_1_project.html", "keys": ["rodin::geometry::polytope::project", "geometry::polytope::project", "polytope::project", "project"]}, {"name": "Rodin::Geometry::Polytope::Traits", "url": "struct_rodin_1_1_geometry_1_1_polytope_1_1_traits.html", "keys": ["rodin::geometry::polytope::traits", "geometry::polytope::traits", "polytope::traits", "traits"]}, {"name": "Rodin::Geometry::Polytope::Traits::HalfSpace", "url": "struct_rodin_1_1_geometry_1_1_polytope_1_1_traits_1_1_half_space.html", "keys": ["rodin::geometry::polytope::traits::halfspace", "geometry::polytope::traits::halfspace", "polytope::traits::halfspace", "traits::halfspace", "halfspace"]}, {"name": "Rodin::Geometry::Shard::PolytopeMap", "url": "struct_rodin_1_1_geometry_1_1_shard_1_1_polytope_map.html", "keys": ["rodin::geometry::shard::polytopemap", "geometry::shard::polytopemap", "shard::polytopemap", "polytopemap"]}, {"name": "Rodin::Geometry::SubMeshBase::PolytopeMap", "url": "struct_rodin_1_1_geometry_1_1_sub_mesh_base_1_1_polytope_map.html", "keys": ["rodin::geometry::submeshbase::polytopemap", "geometry::submeshbase::polytopemap", "submeshbase::polytopemap", "polytopemap"]}, {"name": "Rodin::Heart::CCMLC2014::Model::EvalDataT", "url": "struct_rodin_1_1_heart_1_1_c_c_m_l_c2014_1_1_model_1_1_eval_data_t.html", "keys": ["rodin::heart::ccmlc2014::model::evaldatat", "heart::ccmlc2014::model::evaldatat", "ccmlc2014::model::evaldatat", "model::evaldatat", "evaldatat"]}, {"name": "Rodin::Heart::CCMLC2014::Model::HistoryT", "url": "struct_rodin_1_1_heart_1_1_c_c_m_l_c2014_1_1_model_1_1_history_t.html", "keys": ["rodin::heart::ccmlc2014::model::historyt", "heart::ccmlc2014::model::historyt", "ccmlc2014::model::historyt", "model::historyt", "historyt"]}, {"name": "Rodin::Heart::CCMLC2014::Model::InputT", "url": "struct_rodin_1_1_heart_1_1_c_c_m_l_c2014_1_1_model_1_1_input_t.html", "keys": ["rodin::heart::ccmlc2014::model::inputt", "heart::ccmlc2014::model::inputt", "ccmlc2014::model::inputt", "model::inputt", "inputt"]}, {"name": "Rodin::Heart::CCMLC2014::Model::LocalActiveDataT", "url": "struct_rodin_1_1_heart_1_1_c_c_m_l_c2014_1_1_model_1_1_local_active_data_t.html", "keys": ["rodin::heart::ccmlc2014::model::localactivedatat", "heart::ccmlc2014::model::localactivedatat", "ccmlc2014::model::localactivedatat", "model::localactivedatat", "localactivedatat"]}, {"name": "Rodin::Heart::CCMLC2014::Model::ReportT", "url": "struct_rodin_1_1_heart_1_1_c_c_m_l_c2014_1_1_model_1_1_report_t.html", "keys": ["rodin::heart::ccmlc2014::model::reportt", "heart::ccmlc2014::model::reportt", "ccmlc2014::model::reportt", "model::reportt", "reportt"]}, {"name": "Rodin::Heart::CCMLC2014::Model::StateT", "url": "struct_rodin_1_1_heart_1_1_c_c_m_l_c2014_1_1_model_1_1_state_t.html", "keys": ["rodin::heart::ccmlc2014::model::statet", "heart::ccmlc2014::model::statet", "ccmlc2014::model::statet", "model::statet", "statet"]}, {"name": "Rodin::Heart::CCMLC2014::Physics::Rheology::CarreauYasuda", "url": "struct_rodin_1_1_heart_1_1_c_c_m_l_c2014_1_1_physics_1_1_rheology_1_1_carreau_yasuda.html", "keys": ["rodin::heart::ccmlc2014::physics::rheology::carreauyasuda", "heart::ccmlc2014::physics::rheology::carreauyasuda", "ccmlc2014::physics::rheology::carreauyasuda", "physics::rheology::carreauyasuda", "rheology::carreauyasuda", "carreauyasuda"]}, {"name": "Rodin::Heart::CCMLC2014::Physics::Rheology::Cross", "url": "struct_rodin_1_1_heart_1_1_c_c_m_l_c2014_1_1_physics_1_1_rheology_1_1_cross.html", "keys": ["rodin::heart::ccmlc2014::physics::rheology::cross", "heart::ccmlc2014::physics::rheology::cross", "ccmlc2014::physics::rheology::cross", "physics::rheology::cross", "rheology::cross", "cross"]}, {"name": "Rodin::Heart::CCMLC2014::Physics::Rheology::Newtonian", "url": "struct_rodin_1_1_heart_1_1_c_c_m_l_c2014_1_1_physics_1_1_rheology_1_1_newtonian.html", "keys": ["rodin::heart::ccmlc2014::physics::rheology::newtonian", "heart::ccmlc2014::physics::rheology::newtonian", "ccmlc2014::physics::rheology::newtonian", "physics::rheology::newtonian", "rheology::newtonian", "newtonian"]}, {"name": "Rodin::Heart::CCMLC2014::Physics::Rheology::PowerLaw", "url": "struct_rodin_1_1_heart_1_1_c_c_m_l_c2014_1_1_physics_1_1_rheology_1_1_power_law.html", "keys": ["rodin::heart::ccmlc2014::physics::rheology::powerlaw", "heart::ccmlc2014::physics::rheology::powerlaw", "ccmlc2014::physics::rheology::powerlaw", "physics::rheology::powerlaw", "rheology::powerlaw", "powerlaw"]}, {"name": "Rodin::Heart::CCMLC2014::Physics::Rheology::Quemada", "url": "struct_rodin_1_1_heart_1_1_c_c_m_l_c2014_1_1_physics_1_1_rheology_1_1_quemada.html", "keys": ["rodin::heart::ccmlc2014::physics::rheology::quemada", "heart::ccmlc2014::physics::rheology::quemada", "ccmlc2014::physics::rheology::quemada", "physics::rheology::quemada", "rheology::quemada", "quemada"]}, {"name": "Rodin::Heart::CCMLC2014::Physics::SolverConfig", "url": "struct_rodin_1_1_heart_1_1_c_c_m_l_c2014_1_1_physics_1_1_solver_config.html", "keys": ["rodin::heart::ccmlc2014::physics::solverconfig", "heart::ccmlc2014::physics::solverconfig", "ccmlc2014::physics::solverconfig", "physics::solverconfig", "solverconfig"]}, {"name": "Rodin::Heart::CCMLC2014PassiveLaw", "url": "struct_rodin_1_1_heart_1_1_c_c_m_l_c2014_passive_law.html", "keys": ["rodin::heart::ccmlc2014passivelaw", "heart::ccmlc2014passivelaw", "ccmlc2014passivelaw"]}, {"name": "Rodin::Heart::HolzapfelReducedLaw::Parameters", "url": "struct_rodin_1_1_heart_1_1_holzapfel_reduced_law_1_1_parameters.html", "keys": ["rodin::heart::holzapfelreducedlaw::parameters", "heart::holzapfelreducedlaw::parameters", "holzapfelreducedlaw::parameters", "parameters"]}, {"name": "Rodin::Heart::PassiveEnergyDerivatives", "url": "struct_rodin_1_1_heart_1_1_passive_energy_derivatives.html", "keys": ["rodin::heart::passiveenergyderivatives", "heart::passiveenergyderivatives", "passiveenergyderivatives"]}, {"name": "Rodin::Heart::ReducedInvariantGradient", "url": "struct_rodin_1_1_heart_1_1_reduced_invariant_gradient.html", "keys": ["rodin::heart::reducedinvariantgradient", "heart::reducedinvariantgradient", "reducedinvariantgradient"]}, {"name": "Rodin::Heart::ReducedInvariantHessian", "url": "struct_rodin_1_1_heart_1_1_reduced_invariant_hessian.html", "keys": ["rodin::heart::reducedinvarianthessian", "heart::reducedinvarianthessian", "reducedinvarianthessian"]}, {"name": "Rodin::Heart::ReducedInvariants", "url": "struct_rodin_1_1_heart_1_1_reduced_invariants.html", "keys": ["rodin::heart::reducedinvariants", "heart::reducedinvariants", "reducedinvariants"]}, {"name": "Rodin::IO::HDF5::TransformationData", "url": "struct_rodin_1_1_i_o_1_1_h_d_f5_1_1_transformation_data.html", "keys": ["rodin::io::hdf5::transformationdata", "io::hdf5::transformationdata", "hdf5::transformationdata", "transformationdata"]}, {"name": "Rodin::IO::HDF5::XDMFTopologyLayout", "url": "struct_rodin_1_1_i_o_1_1_h_d_f5_1_1_x_d_m_f_topology_layout.html", "keys": ["rodin::io::hdf5::xdmftopologylayout", "io::hdf5::xdmftopologylayout", "hdf5::xdmftopologylayout", "xdmftopologylayout"]}, {"name": "Rodin::IO::MEDIT::ParseEntity::Data", "url": "struct_rodin_1_1_i_o_1_1_m_e_d_i_t_1_1_parse_entity_1_1_data.html", "keys": ["rodin::io::medit::parseentity::data", "io::medit::parseentity::data", "medit::parseentity::data", "parseentity::data", "data"]}, {"name": "Rodin::IO::MFEM::FaceOrderEntry", "url": "struct_rodin_1_1_i_o_1_1_m_f_e_m_1_1_face_order_entry.html", "keys": ["rodin::io::mfem::faceorderentry", "io::mfem::faceorderentry", "mfem::faceorderentry", "faceorderentry"]}, {"name": "Rodin::IO::MFEM::GridFunctionHeader", "url": "struct_rodin_1_1_i_o_1_1_m_f_e_m_1_1_grid_function_header.html", "keys": ["rodin::io::mfem::gridfunctionheader", "io::mfem::gridfunctionheader", "mfem::gridfunctionheader", "gridfunctionheader"]}, {"name": "Rodin::IO::MFEM::MeshHeader", "url": "struct_rodin_1_1_i_o_1_1_m_f_e_m_1_1_mesh_header.html", "keys": ["rodin::io::mfem::meshheader", "io::mfem::meshheader", "mfem::meshheader", "meshheader"]}, {"name": "Rodin::IO::MFEM::MeshVersion", "url": "struct_rodin_1_1_i_o_1_1_m_f_e_m_1_1_mesh_version.html", "keys": ["rodin::io::mfem::meshversion", "io::mfem::meshversion", "mfem::meshversion", "meshversion"]}, {"name": "Rodin::IO::MFEM::ParseEmptyLine", "url": "struct_rodin_1_1_i_o_1_1_m_f_e_m_1_1_parse_empty_line.html", "keys": ["rodin::io::mfem::parseemptyline", "io::mfem::parseemptyline", "mfem::parseemptyline", "parseemptyline"]}, {"name": "Rodin::IO::MFEM::ParseEmptyLineOrComment", "url": "struct_rodin_1_1_i_o_1_1_m_f_e_m_1_1_parse_empty_line_or_comment.html", "keys": ["rodin::io::mfem::parseemptylineorcomment", "io::mfem::parseemptylineorcomment", "mfem::parseemptylineorcomment", "parseemptylineorcomment"]}, {"name": "Rodin::IO::MFEM::ParseGeometry", "url": "struct_rodin_1_1_i_o_1_1_m_f_e_m_1_1_parse_geometry.html", "keys": ["rodin::io::mfem::parsegeometry", "io::mfem::parsegeometry", "mfem::parsegeometry", "parsegeometry"]}, {"name": "Rodin::IO::MFEM::ParseGeometry::Data", "url": "struct_rodin_1_1_i_o_1_1_m_f_e_m_1_1_parse_geometry_1_1_data.html", "keys": ["rodin::io::mfem::parsegeometry::data", "io::mfem::parsegeometry::data", "mfem::parsegeometry::data", "parsegeometry::data", "data"]}, {"name": "Rodin::IO::MFEM::ParseKeyword", "url": "struct_rodin_1_1_i_o_1_1_m_f_e_m_1_1_parse_keyword.html", "keys": ["rodin::io::mfem::parsekeyword", "io::mfem::parsekeyword", "mfem::parsekeyword", "parsekeyword"]}, {"name": "Rodin::IO::MFEM::ParseUnsignedInteger", "url": "struct_rodin_1_1_i_o_1_1_m_f_e_m_1_1_parse_unsigned_integer.html", "keys": ["rodin::io::mfem::parseunsignedinteger", "io::mfem::parseunsignedinteger", "mfem::parseunsignedinteger", "parseunsignedinteger"]}, {"name": "Rodin::IO::XDMF::FilePatterns", "url": "struct_rodin_1_1_i_o_1_1_x_d_m_f_1_1_file_patterns.html", "keys": ["rodin::io::xdmf::filepatterns", "io::xdmf::filepatterns", "xdmf::filepatterns", "filepatterns"]}, {"name": "Rodin::IO::XDMF::GridOptions", "url": "struct_rodin_1_1_i_o_1_1_x_d_m_f_1_1_grid_options.html", "keys": ["rodin::io::xdmf::gridoptions", "io::xdmf::gridoptions", "xdmf::gridoptions", "gridoptions"]}, {"name": "Rodin::IO::XDMF::Keyword", "url": "struct_rodin_1_1_i_o_1_1_x_d_m_f_1_1_keyword.html", "keys": ["rodin::io::xdmf::keyword", "io::xdmf::keyword", "xdmf::keyword", "keyword"]}, {"name": "Rodin::IO::XDMF::SnapshotRecord::AttributeRecord", "url": "struct_rodin_1_1_i_o_1_1_x_d_m_f_1_1_snapshot_record_1_1_attribute_record.html", "keys": ["rodin::io::xdmf::snapshotrecord::attributerecord", "io::xdmf::snapshotrecord::attributerecord", "xdmf::snapshotrecord::attributerecord", "snapshotrecord::attributerecord", "attributerecord"]}, {"name": "Rodin::IO::XDMF::SnapshotRecord::MeshAttributeRecord", "url": "struct_rodin_1_1_i_o_1_1_x_d_m_f_1_1_snapshot_record_1_1_mesh_attribute_record.html", "keys": ["rodin::io::xdmf::snapshotrecord::meshattributerecord", "io::xdmf::snapshotrecord::meshattributerecord", "xdmf::snapshotrecord::meshattributerecord", "snapshotrecord::meshattributerecord", "meshattributerecord"]}, {"name": "Rodin::IO::XDMF::SnapshotRecord::PieceMeta", "url": "struct_rodin_1_1_i_o_1_1_x_d_m_f_1_1_snapshot_record_1_1_piece_meta.html", "keys": ["rodin::io::xdmf::snapshotrecord::piecemeta", "io::xdmf::snapshotrecord::piecemeta", "xdmf::snapshotrecord::piecemeta", "snapshotrecord::piecemeta", "piecemeta"]}, {"name": "Rodin::IndexArrayCompare", "url": "struct_rodin_1_1_index_array_compare.html", "keys": ["rodin::indexarraycompare", "indexarraycompare"]}, {"name": "Rodin::IndexArrayEquality", "url": "struct_rodin_1_1_index_array_equality.html", "keys": ["rodin::indexarrayequality", "indexarrayequality"]}, {"name": "Rodin::IndexSetEquality", "url": "struct_rodin_1_1_index_set_equality.html", "keys": ["rodin::indexsetequality", "indexsetequality"]}, {"name": "Rodin::IndexSetHash", "url": "struct_rodin_1_1_index_set_hash.html", "keys": ["rodin::indexsethash", "indexsethash"]}, {"name": "Rodin::MMG::NoSplitT", "url": "struct_rodin_1_1_m_m_g_1_1_no_split_t.html", "keys": ["rodin::mmg::nosplitt", "mmg::nosplitt", "nosplitt"]}, {"name": "Rodin::MMG::Split", "url": "struct_rodin_1_1_m_m_g_1_1_split.html", "keys": ["rodin::mmg::split", "mmg::split", "split"]}, {"name": "Rodin::Math::RungeKutta::RK2", "url": "struct_rodin_1_1_math_1_1_runge_kutta_1_1_r_k2.html", "keys": ["rodin::math::rungekutta::rk2", "math::rungekutta::rk2", "rungekutta::rk2", "rk2"]}, {"name": "Rodin::Math::RungeKutta::RK4", "url": "struct_rodin_1_1_math_1_1_runge_kutta_1_1_r_k4.html", "keys": ["rodin::math::rungekutta::rk4", "math::rungekutta::rk4", "rungekutta::rk4", "rk4"]}, {"name": "Rodin::PETSc::Assembly::MatrixSetup::Options", "url": "struct_rodin_1_1_p_e_t_sc_1_1_assembly_1_1_matrix_setup_1_1_options.html", "keys": ["rodin::petsc::assembly::matrixsetup::options", "petsc::assembly::matrixsetup::options", "assembly::matrixsetup::options", "matrixsetup::options", "options"]}, {"name": "Rodin::PETSc::Assembly::VectorSetup::Options", "url": "struct_rodin_1_1_p_e_t_sc_1_1_assembly_1_1_vector_setup_1_1_options.html", "keys": ["rodin::petsc::assembly::vectorsetup::options", "petsc::assembly::vectorsetup::options", "assembly::vectorsetup::options", "vectorsetup::options", "options"]}, {"name": "Rodin::QF::QuadratureFormula::Lookup::Entry", "url": "struct_rodin_1_1_q_f_1_1_quadrature_formula_1_1_lookup_1_1_entry.html", "keys": ["rodin::qf::quadratureformula::lookup::entry", "qf::quadratureformula::lookup::entry", "quadratureformula::lookup::entry", "lookup::entry", "entry"]}, {"name": "Rodin::QF::QuadratureFormulaBase::Samples", "url": "struct_rodin_1_1_q_f_1_1_quadrature_formula_base_1_1_samples.html", "keys": ["rodin::qf::quadratureformulabase::samples", "qf::quadratureformulabase::samples", "quadratureformulabase::samples", "samples"]}, {"name": "Rodin::Solid::ActiveFiberLaw::Parameters", "url": "struct_rodin_1_1_solid_1_1_active_fiber_law_1_1_parameters.html", "keys": ["rodin::solid::activefiberlaw::parameters", "solid::activefiberlaw::parameters", "activefiberlaw::parameters", "parameters"]}, {"name": "Rodin::Solid::ActiveFiberLaw::Parameters::Initial", "url": "struct_rodin_1_1_solid_1_1_active_fiber_law_1_1_parameters_1_1_initial.html", "keys": ["rodin::solid::activefiberlaw::parameters::initial", "solid::activefiberlaw::parameters::initial", "activefiberlaw::parameters::initial", "parameters::initial", "initial"]}, {"name": "Rodin::Solid::ActiveFiberLaw::Residual", "url": "struct_rodin_1_1_solid_1_1_active_fiber_law_1_1_residual.html", "keys": ["rodin::solid::activefiberlaw::residual", "solid::activefiberlaw::residual", "activefiberlaw::residual", "residual"]}, {"name": "Rodin::Solid::ActiveFiberLaw::Response", "url": "struct_rodin_1_1_solid_1_1_active_fiber_law_1_1_response.html", "keys": ["rodin::solid::activefiberlaw::response", "solid::activefiberlaw::response", "activefiberlaw::response", "response"]}, {"name": "Rodin::Solid::ActiveFiberLaw::State", "url": "struct_rodin_1_1_solid_1_1_active_fiber_law_1_1_state.html", "keys": ["rodin::solid::activefiberlaw::state", "solid::activefiberlaw::state", "activefiberlaw::state", "state"]}, {"name": "Rodin::Solid::HolzapfelOgden::Parameters", "url": "struct_rodin_1_1_solid_1_1_holzapfel_ogden_1_1_parameters.html", "keys": ["rodin::solid::holzapfelogden::parameters", "solid::holzapfelogden::parameters", "holzapfelogden::parameters", "parameters"]}, {"name": "Rodin::Solid::Internal::FaceKinematics", "url": "struct_rodin_1_1_solid_1_1_internal_1_1_face_kinematics.html", "keys": ["rodin::solid::internal::facekinematics", "solid::internal::facekinematics", "internal::facekinematics", "facekinematics"]}, {"name": "Rodin::Solid::Tags::Activation", "url": "struct_rodin_1_1_solid_1_1_tags_1_1_activation.html", "keys": ["rodin::solid::tags::activation", "solid::tags::activation", "tags::activation", "activation"]}, {"name": "Rodin::Solid::Tags::ActiveBeta", "url": "struct_rodin_1_1_solid_1_1_tags_1_1_active_beta.html", "keys": ["rodin::solid::tags::activebeta", "solid::tags::activebeta", "tags::activebeta", "activebeta"]}, {"name": "Rodin::Solid::Tags::ActiveExtension", "url": "struct_rodin_1_1_solid_1_1_tags_1_1_active_extension.html", "keys": ["rodin::solid::tags::activeextension", "solid::tags::activeextension", "tags::activeextension", "activeextension"]}, {"name": "Rodin::Solid::Tags::ActiveGamma", "url": "struct_rodin_1_1_solid_1_1_tags_1_1_active_gamma.html", "keys": ["rodin::solid::tags::activegamma", "solid::tags::activegamma", "tags::activegamma", "activegamma"]}, {"name": "Rodin::Solid::Tags::CellIndex", "url": "struct_rodin_1_1_solid_1_1_tags_1_1_cell_index.html", "keys": ["rodin::solid::tags::cellindex", "solid::tags::cellindex", "tags::cellindex", "cellindex"]}, {"name": "Rodin::Solid::Tags::ElectricalActivation", "url": "struct_rodin_1_1_solid_1_1_tags_1_1_electrical_activation.html", "keys": ["rodin::solid::tags::electricalactivation", "solid::tags::electricalactivation", "tags::electricalactivation", "electricalactivation"]}, {"name": "Rodin::Solid::Tags::FiberDirection", "url": "struct_rodin_1_1_solid_1_1_tags_1_1_fiber_direction.html", "keys": ["rodin::solid::tags::fiberdirection", "solid::tags::fiberdirection", "tags::fiberdirection", "fiberdirection"]}, {"name": "Rodin::Solid::Tags::PreviousActiveBeta", "url": "struct_rodin_1_1_solid_1_1_tags_1_1_previous_active_beta.html", "keys": ["rodin::solid::tags::previousactivebeta", "solid::tags::previousactivebeta", "tags::previousactivebeta", "previousactivebeta"]}, {"name": "Rodin::Solid::Tags::PreviousActiveExtension", "url": "struct_rodin_1_1_solid_1_1_tags_1_1_previous_active_extension.html", "keys": ["rodin::solid::tags::previousactiveextension", "solid::tags::previousactiveextension", "tags::previousactiveextension", "previousactiveextension"]}, {"name": "Rodin::Solid::Tags::PreviousActiveGamma", "url": "struct_rodin_1_1_solid_1_1_tags_1_1_previous_active_gamma.html", "keys": ["rodin::solid::tags::previousactivegamma", "solid::tags::previousactivegamma", "tags::previousactivegamma", "previousactivegamma"]}, {"name": "Rodin::Solid::Tags::PreviousFiberStrain", "url": "struct_rodin_1_1_solid_1_1_tags_1_1_previous_fiber_strain.html", "keys": ["rodin::solid::tags::previousfiberstrain", "solid::tags::previousfiberstrain", "tags::previousfiberstrain", "previousfiberstrain"]}, {"name": "Rodin::Solid::Tags::QuadraturePointIndex", "url": "struct_rodin_1_1_solid_1_1_tags_1_1_quadrature_point_index.html", "keys": ["rodin::solid::tags::quadraturepointindex", "solid::tags::quadraturepointindex", "tags::quadraturepointindex", "quadraturepointindex"]}, {"name": "Rodin::Solid::Tags::SheetDirection", "url": "struct_rodin_1_1_solid_1_1_tags_1_1_sheet_direction.html", "keys": ["rodin::solid::tags::sheetdirection", "solid::tags::sheetdirection", "tags::sheetdirection", "sheetdirection"]}, {"name": "Rodin::Solid::Tags::SheetNormalDirection", "url": "struct_rodin_1_1_solid_1_1_tags_1_1_sheet_normal_direction.html", "keys": ["rodin::solid::tags::sheetnormaldirection", "solid::tags::sheetnormaldirection", "tags::sheetnormaldirection", "sheetnormaldirection"]}, {"name": "Rodin::Solid::Tags::TimeStep", "url": "struct_rodin_1_1_solid_1_1_tags_1_1_time_step.html", "keys": ["rodin::solid::tags::timestep", "solid::tags::timestep", "tags::timestep", "timestep"]}, {"name": "Rodin::Solver::Info", "url": "struct_rodin_1_1_solver_1_1_info.html", "keys": ["rodin::solver::info", "solver::info", "info"]}, {"name": "Rodin::Solver::NewtonSolver::Report", "url": "struct_rodin_1_1_solver_1_1_newton_solver_1_1_report.html", "keys": ["rodin::solver::newtonsolver::report", "solver::newtonsolver::report", "newtonsolver::report", "report"]}, {"name": "Rodin::Solver::NewtonSolver::StepResult", "url": "struct_rodin_1_1_solver_1_1_newton_solver_1_1_step_result.html", "keys": ["rodin::solver::newtonsolver::stepresult", "solver::newtonsolver::stepresult", "newtonsolver::stepresult", "stepresult"]}, {"name": "Rodin::Test::FDProbeReport", "url": "struct_rodin_1_1_test_1_1_f_d_probe_report.html", "keys": ["rodin::test::fdprobereport", "test::fdprobereport", "fdprobereport"]}, {"name": "Rodin::Test::Random::detail::UniformIntDistType", "url": "struct_rodin_1_1_test_1_1_random_1_1detail_1_1_uniform_int_dist_type.html", "keys": ["rodin::test::random::detail::uniformintdisttype", "test::random::detail::uniformintdisttype", "random::detail::uniformintdisttype", "detail::uniformintdisttype", "uniformintdisttype"]}, {"name": "Rodin::Utility::BottomTemplate", "url": "struct_rodin_1_1_utility_1_1_bottom_template.html", "keys": ["rodin::utility::bottomtemplate", "utility::bottomtemplate", "bottomtemplate"]}, {"name": "Rodin::Utility::DependentValue", "url": "struct_rodin_1_1_utility_1_1_dependent_value.html", "keys": ["rodin::utility::dependentvalue", "utility::dependentvalue", "dependentvalue"]}, {"name": "Rodin::Utility::Extract", "url": "struct_rodin_1_1_utility_1_1_extract.html", "keys": ["rodin::utility::extract", "utility::extract", "extract"]}, {"name": "Rodin::Utility::HasTypeMember", "url": "struct_rodin_1_1_utility_1_1_has_type_member.html", "keys": ["rodin::utility::hastypemember", "utility::hastypemember", "hastypemember"]}, {"name": "Rodin::Utility::HasValueMember", "url": "struct_rodin_1_1_utility_1_1_has_value_member.html", "keys": ["rodin::utility::hasvaluemember", "utility::hasvaluemember", "hasvaluemember"]}, {"name": "Rodin::Utility::Internal::AllImpl", "url": "struct_rodin_1_1_utility_1_1_internal_1_1_all_impl.html", "keys": ["rodin::utility::internal::allimpl", "utility::internal::allimpl", "internal::allimpl", "allimpl"]}, {"name": "Rodin::Utility::Internal::AtImpl", "url": "struct_rodin_1_1_utility_1_1_internal_1_1_at_impl.html", "keys": ["rodin::utility::internal::atimpl", "utility::internal::atimpl", "internal::atimpl", "atimpl"]}, {"name": "Rodin::Utility::IsOneOf", "url": "struct_rodin_1_1_utility_1_1_is_one_of.html", "keys": ["rodin::utility::isoneof", "utility::isoneof", "isoneof"]}, {"name": "Rodin::Utility::IsSpecialization", "url": "struct_rodin_1_1_utility_1_1_is_specialization.html", "keys": ["rodin::utility::isspecialization", "utility::isspecialization", "isspecialization"]}, {"name": "Rodin::Utility::Make", "url": "struct_rodin_1_1_utility_1_1_make.html", "keys": ["rodin::utility::make", "utility::make", "make"]}, {"name": "Rodin::Utility::Overloaded", "url": "struct_rodin_1_1_utility_1_1_overloaded.html", "keys": ["rodin::utility::overloaded", "utility::overloaded", "overloaded"]}, {"name": "Rodin::Utility::Product", "url": "struct_rodin_1_1_utility_1_1_product.html", "keys": ["rodin::utility::product", "utility::product", "product"]}, {"name": "Rodin::Utility::Repeat", "url": "struct_rodin_1_1_utility_1_1_repeat.html", "keys": ["rodin::utility::repeat", "utility::repeat", "repeat"]}, {"name": "Rodin::Utility::UnwrapRefDecay", "url": "struct_rodin_1_1_utility_1_1_unwrap_ref_decay.html", "keys": ["rodin::utility::unwraprefdecay", "utility::unwraprefdecay", "unwraprefdecay"]}, {"name": "Rodin::Utility::UnwrapReference", "url": "struct_rodin_1_1_utility_1_1_unwrap_reference.html", "keys": ["rodin::utility::unwrapreference", "utility::unwrapreference", "unwrapreference"]}, {"name": "Rodin::Utility::Zip", "url": "struct_rodin_1_1_utility_1_1_zip.html", "keys": ["rodin::utility::zip", "utility::zip", "zip"]}, {"name": "Rodin::Variational::AllPETScTrialOrTest", "url": "struct_rodin_1_1_variational_1_1_all_p_e_t_sc_trial_or_test.html", "keys": ["rodin::variational::allpetsctrialortest", "variational::allpetsctrialortest", "allpetsctrialortest"]}, {"name": "Rodin::Variational::BoundaryHit", "url": "struct_rodin_1_1_variational_1_1_boundary_hit.html", "keys": ["rodin::variational::boundaryhit", "variational::boundaryhit", "boundaryhit"]}, {"name": "Rodin::Variational::DubinerTetrahedron::Gradients", "url": "struct_rodin_1_1_variational_1_1_dubiner_tetrahedron_1_1_gradients.html", "keys": ["rodin::variational::dubinertetrahedron::gradients", "variational::dubinertetrahedron::gradients", "dubinertetrahedron::gradients", "gradients"]}, {"name": "Rodin::Variational::DubinerTetrahedron::Values", "url": "struct_rodin_1_1_variational_1_1_dubiner_tetrahedron_1_1_values.html", "keys": ["rodin::variational::dubinertetrahedron::values", "variational::dubinertetrahedron::values", "dubinertetrahedron::values", "values"]}, {"name": "Rodin::Variational::FiniteElementSpaceBase::LocalIndex", "url": "struct_rodin_1_1_variational_1_1_finite_element_space_base_1_1_local_index.html", "keys": ["rodin::variational::finiteelementspacebase::localindex", "variational::finiteelementspacebase::localindex", "finiteelementspacebase::localindex", "localindex"]}, {"name": "Rodin::Variational::H1MassOperands", "url": "struct_rodin_1_1_variational_1_1_h1_mass_operands.html", "keys": ["rodin::variational::h1massoperands", "variational::h1massoperands", "h1massoperands"]}, {"name": "Rodin::Variational::IsPETScTestFunction", "url": "struct_rodin_1_1_variational_1_1_is_p_e_t_sc_test_function.html", "keys": ["rodin::variational::ispetsctestfunction", "variational::ispetsctestfunction", "ispetsctestfunction"]}, {"name": "Rodin::Variational::IsPETScTrialFunction", "url": "struct_rodin_1_1_variational_1_1_is_p_e_t_sc_trial_function.html", "keys": ["rodin::variational::ispetsctrialfunction", "variational::ispetsctrialfunction", "ispetsctrialfunction"]}, {"name": "Rodin::Variational::IsTestFunction", "url": "struct_rodin_1_1_variational_1_1_is_test_function.html", "keys": ["rodin::variational::istestfunction", "variational::istestfunction", "istestfunction"]}, {"name": "Rodin::Variational::IsTrialFunction", "url": "struct_rodin_1_1_variational_1_1_is_trial_function.html", "keys": ["rodin::variational::istrialfunction", "variational::istrialfunction", "istrialfunction"]}, {"name": "Rodin::Variational::PyramidIndex", "url": "struct_rodin_1_1_variational_1_1_pyramid_index.html", "keys": ["rodin::variational::pyramidindex", "variational::pyramidindex", "pyramidindex"]}, {"name": "Rodin::Variational::PyramidIndex::IJ", "url": "struct_rodin_1_1_variational_1_1_pyramid_index_1_1_i_j.html", "keys": ["rodin::variational::pyramidindex::ij", "variational::pyramidindex::ij", "pyramidindex::ij", "ij"]}, {"name": "Rodin::Variational::PyramidModal::Gradients", "url": "struct_rodin_1_1_variational_1_1_pyramid_modal_1_1_gradients.html", "keys": ["rodin::variational::pyramidmodal::gradients", "variational::pyramidmodal::gradients", "pyramidmodal::gradients", "gradients"]}, {"name": "Rodin::Variational::PyramidModal::Values", "url": "struct_rodin_1_1_variational_1_1_pyramid_modal_1_1_values.html", "keys": ["rodin::variational::pyramidmodal::values", "variational::pyramidmodal::values", "pyramidmodal::values", "values"]}, {"name": "Rodin::Variational::ReferenceTabulation::Gradients", "url": "struct_rodin_1_1_variational_1_1_reference_tabulation_1_1_gradients.html", "keys": ["rodin::variational::referencetabulation::gradients", "variational::referencetabulation::gradients", "referencetabulation::gradients", "gradients"]}, {"name": "Rodin::Variational::ReferenceTabulation::Values", "url": "struct_rodin_1_1_variational_1_1_reference_tabulation_1_1_values.html", "keys": ["rodin::variational::referencetabulation::values", "variational::referencetabulation::values", "referencetabulation::values", "values"]}]);
