/* In-page navigation for the generated documentation. */
(() => {
  const content = document.querySelector('main > article > .m-container > .m-row > .m-col-l-10');
  if (!content) return;

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

  const navigation = document.createElement('nav');
  navigation.className = 'rodin-page-nav';
  navigation.setAttribute('aria-label', 'On this page');
  const title = document.createElement('h2');
  title.textContent = 'On this page';
  navigation.append(title);
  const list = document.createElement('ol');
  const parents = [{ level: 0, list }];
  for (const entry of entries) {
    const level = Number(entry.heading.tagName.slice(1));
    while (parents[parents.length - 1].level >= level) parents.pop();
    const parent = parents[parents.length - 1];
    if (!parent.list) {
      parent.list = document.createElement('ol');
      parent.item.append(parent.list);
    }
    const item = document.createElement('li');
    entry.link = document.createElement('a');
    entry.link.href = '#' + encodeURIComponent(entry.target);
    entry.link.textContent = entry.title;
    item.append(entry.link);
    parent.list.append(item);
    parents.push({ level, item });
  }
  navigation.append(list);
  content.parentElement.prepend(navigation);
  const layout = () => {
    const unit = parseFloat(getComputedStyle(document.documentElement).fontSize);
    const gutter = content.getBoundingClientRect().left;
    const width = Math.max(0, Math.min(13 * unit, gutter - 2 * unit));
    navigation.style.width = width + 'px';
    navigation.style.left = (gutter - width - unit) + 'px';
    navigation.classList.toggle('rodin-page-nav-fit', width >= 9 * unit);
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
    if (Math.ceil(window.scrollY + window.innerHeight) >= document.documentElement.scrollHeight) {
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
  const relayout = () => {
    layout();
    schedule();
  };
  window.addEventListener('resize', relayout, { passive: true });
  window.addEventListener('load', relayout, { once: true });
  layout();
  update();
})();
