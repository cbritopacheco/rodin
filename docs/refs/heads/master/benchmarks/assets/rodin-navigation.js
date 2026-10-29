/* In-page navigation for the generated documentation. */
(() => {
  const benchmark = document.querySelector('[data-rodin-content]');
  const content = benchmark || document.querySelector('main > article > .m-container > .m-row > .m-col-l-10');
  if (!content) return;

  // The homepage signature already identifies the project.
  if (content.querySelector(':scope > .rodin-main-brand')) {
    content.querySelector(':scope > h1')?.remove();
  }

  const prepareContents = toc => {
    for (const item of toc.querySelectorAll('li')) {
      const children = item.querySelector(':scope > ul, :scope > ol');
      const link = item.querySelector(':scope > a');
      if (!children || !link || item.querySelector(':scope > button')) continue;
      const toggle = document.createElement('button');
      toggle.type = 'button';
      toggle.className = 'rodin-toc-toggle';
      toggle.textContent = '−';
      toggle.setAttribute('aria-expanded', 'true');
      toggle.setAttribute('aria-label', 'Collapse ' + link.textContent.trim());
      item.classList.add('rodin-toc-branch');
      item.insertBefore(toggle, link);
      toggle.addEventListener('click', () => {
        children.hidden = !children.hidden;
        toggle.textContent = children.hidden ? '+' : '−';
        toggle.setAttribute('aria-expanded', String(!children.hidden));
        toggle.setAttribute('aria-label',
          (children.hidden ? 'Expand ' : 'Collapse ') + link.textContent.trim());
      });
    }
  };
  for (const toc of content.querySelectorAll('nav, .m-table-of-contents')) {
    const heading = toc.querySelector(':scope > h2, :scope > h3');
    if (!toc.matches('[aria-label*="table of contents"], .m-table-of-contents')
        && heading?.textContent.trim() !== 'Contents') continue;
    toc.classList.add('rodin-toc');
    const disclosure = document.createElement('details');
    disclosure.open = true;
    const summary = document.createElement('summary');
    summary.textContent = 'Contents';
    heading?.remove();
    disclosure.append(summary, ...toc.childNodes);
    toc.append(disclosure);
    prepareContents(toc);
  }

  let selected = null;
  const entries = [];
  const targets = new Set();
  for (const heading of content.querySelectorAll('h1, h2, h3')) {
    if (heading.closest('nav, .m-doc-details')) continue;
    const title = heading.tagName === 'H1' ? 'Overview' : heading.textContent.replace(/¶/g, '').trim();
    if (!title) continue;
    const section = heading.closest('section[id]');
    const sectionHeading = section?.querySelector('h1, h2, h3');
    let target = heading.id || (sectionHeading === heading ? section.id : null);
    if (!target) {
      let index = entries.length + 1;
      do { target = `rodin-section-${index++}`; } while (document.getElementById(target));
      heading.id = target;
    }
    if (targets.has(target)) continue;
    targets.add(target);
    entries.push({ heading, target, title });
  }
  if (entries.length < 2) return;

  const positionKey = 'rodin-sidebar-side';
  let side = 'left';
  try { if (localStorage.getItem(positionKey) === 'right') side = 'right'; } catch {}
  const navigation = document.createElement('nav');
  navigation.className = 'rodin-page-nav';
  navigation.setAttribute('aria-label', 'On this page');
  const title = document.createElement('h2');
  title.textContent = 'On this page';
  const header = document.createElement('summary');
  header.className = 'rodin-page-nav-header';
  const positions = document.createElement('div');
  positions.className = 'rodin-page-nav-positions';
  positions.setAttribute('role', 'group');
  positions.setAttribute('aria-label', 'Sidebar position');
  const controls = [];
  for (const location of ['left', 'right']) {
    const button = document.createElement('button');
    button.type = 'button';
    button.dataset.side = location;
    button.setAttribute('aria-label', 'Sidebar on ' + location);
    button.setAttribute('aria-pressed', String(location === side));
    button.title = 'Sidebar on ' + location;
    const icon = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    icon.setAttribute('viewBox', '0 0 24 24');
    icon.setAttribute('aria-hidden', 'true');
    const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    path.setAttribute('d', location === 'left' ? 'M3 4h18v16H3z M9 4v16 M5 8h2 M5 12h2'
      : 'M3 4h18v16H3z M15 4v16 M17 8h2 M17 12h2');
    icon.append(path);
    button.append(icon);
    button.addEventListener('click', () => {
      side = location;
      try { localStorage.setItem(positionKey, side); } catch {}
      for (const control of controls) {
        control.setAttribute('aria-pressed', String(control.dataset.side === side));
      }
      relayout();
    });
    controls.push(button);
    positions.append(button);
  }
  header.append(title, positions);
  const disclosure = document.createElement('details');
  disclosure.open = true;
  disclosure.append(header);
  navigation.append(disclosure);
  const list = document.createElement('ul');
  const parents = [{ level: 0, list }];
  for (const entry of entries) {
    const level = Number(entry.heading.tagName.slice(1));
    while (parents[parents.length - 1].level >= level) parents.pop();
    const parent = parents[parents.length - 1];
    if (!parent.list) {
      parent.list = document.createElement('ul');
      parent.item.append(parent.list);
    }
    const item = document.createElement('li');
    const category = entry.heading.closest('[data-category]')?.dataset.category;
    if (category) item.dataset.category = category;
    entry.link = document.createElement('a');
    entry.link.href = '#' + encodeURIComponent(entry.target);
    entry.link.textContent = entry.title;
    item.append(entry.link);
    parent.list.append(item);
    parents.push({ level, item });
  }
  disclosure.append(list);
  if (benchmark) {
    const contents = document.createElement('nav');
    contents.className = 'm-block m-default';
    contents.setAttribute('aria-label', 'Benchmark table of contents');
    const body = document.createElement('details');
    body.open = true;
    const label = document.createElement('summary');
    label.textContent = 'Contents';
    body.append(label);
    contents.append(body);
    const sections = list.cloneNode(true);
    const overview = sections.firstElementChild;
    if (overview?.querySelector('ul')) {
      body.append(overview.querySelector('ul'));
    } else {
      body.append(sections);
    }
    const configuration = benchmark.querySelector('.benchmark-configuration');
    (configuration || benchmark.querySelector('h1')).after(contents);
    navigation.addEventListener('click', event => {
      const link = event.target.closest('a');
      const entry = entries.find(item => item.link === link);
      if (!entry) return;
      event.preventDefault();
      selected = entry;
      const reduced = matchMedia('(prefers-reduced-motion: reduce)').matches;
      entry.heading.scrollIntoView({ block: 'start', behavior: reduced ? 'instant' : 'smooth' });
      history.replaceState(history.state, '', '#' + encodeURIComponent(entry.target));
      schedule();
    });
  }
  content.parentElement.prepend(navigation);
  prepareContents(navigation);
  for (const toc of content.querySelectorAll('nav[aria-label*="table of contents"]')) {
    prepareContents(toc);
  }
  // The rendered heading is the source of truth, including uncategorized pages.
  const synchronizeColors = () => {
    for (const entry of entries) {
      entry.link.parentElement.style.setProperty('--rodin-section-color',
        getComputedStyle(entry.heading).color);
    }
  };
  synchronizeColors();
  const themeObserver = new MutationObserver(synchronizeColors);
  themeObserver.observe(document.documentElement, {
    attributes: true, attributeFilter: ['data-theme']
  });
  matchMedia('(prefers-color-scheme: light)').addEventListener('change', synchronizeColors);

  const layout = () => {
    if (benchmark) document.body.dataset.contentsSide = side;
    const unit = parseFloat(getComputedStyle(document.documentElement).fontSize);
    const bounds = content.getBoundingClientRect();
    const gutter = side === 'left' ? bounds.left : document.documentElement.clientWidth - bounds.right;
    const width = Math.max(0, Math.min(13 * unit, gutter - 2 * unit));
    const homeWidth = Math.max(0, Math.min(13 * unit, bounds.left - 2 * unit)) + unit;
    document.documentElement.style.setProperty('--rodin-home-width', homeWidth + 'px');
    navigation.style.width = width + 'px';
    navigation.style.left = (side === 'left' ? bounds.left - width - unit : bounds.right + unit) + 'px';
    navigation.dataset.side = side;
    navigation.classList.toggle('rodin-page-nav-fit', width >= 9 * unit);
    const controlGutter = list.getBoundingClientRect().left + .15 * unit;
    for (const toggle of navigation.querySelectorAll('.rodin-toc-toggle')) {
      if (toggle.parentElement.getClientRects().length) {
        toggle.style.left = controlGutter - toggle.parentElement.getBoundingClientRect().left + 'px';
      }
    }
  };

  let current;
  let scheduled = false;
  const update = () => {
    scheduled = false;
    if (!navigation.getClientRects().length) return;
    let active = entries[0];
    for (const entry of entries) {
      if (!entry.heading.getClientRects().length) continue;
      if (entry.heading.getBoundingClientRect().top > 140) break;
      active = entry;
    }
    const bottom = Math.ceil(window.scrollY + window.innerHeight) >= document.documentElement.scrollHeight;
    if (benchmark && selected && bottom) {
      const bounds = selected.heading.getBoundingClientRect();
      // Near the page end a target cannot always reach the reading line.
      if (bounds.top >= 0 && bounds.top < window.innerHeight) active = selected;
    }
    if (!benchmark && Math.ceil(window.scrollY + window.innerHeight) >= document.documentElement.scrollHeight) {
      for (const entry of entries) {
        if (entry.heading.getClientRects().length) active = entry;
      }
    }
    if (current !== active) {
      if (current) current.link.removeAttribute('aria-current');
      active.link.setAttribute('aria-current', 'location');
      current = active;
    }
    // Scroll the navigation list only; never move the document or keyboard focus.
    const link = active.link.getBoundingClientRect();
    const bounds = list.getBoundingClientRect();
    if (link.top < bounds.top) list.scrollTop += link.top - bounds.top;
    else if (link.bottom > bounds.bottom) list.scrollTop += link.bottom - bounds.bottom;
  };
  const schedule = () => {
    if (!scheduled) {
      scheduled = true;
      requestAnimationFrame(update);
    }
  };
  window.addEventListener('scroll', schedule, { passive: true });
  if (benchmark) {
    const resume = () => { selected = null; schedule(); };
    window.addEventListener('wheel', resume, { passive: true });
    window.addEventListener('touchstart', resume, { passive: true });
    window.addEventListener('keydown', event => {
      if (['ArrowUp', 'ArrowDown', 'PageUp', 'PageDown', 'Home', 'End', ' '].includes(event.key)) resume();
    });
    window.addEventListener('pointerdown', event => {
      if (!navigation.contains(event.target)) resume();
    }, { passive: true });
  }
  const relayout = () => {
    layout();
    schedule();
  };
  window.addEventListener('resize', relayout, { passive: true });
  window.addEventListener('load', relayout, { once: true });
  layout();
  update();
})();
