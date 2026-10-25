(() => {
  const key = 'rodin-color-theme';
  const system = window.matchMedia('(prefers-color-scheme: light)');
  const requested = new URLSearchParams(location.search).get('theme');
  let saved = null;
  try { saved = localStorage.getItem(key); } catch {}
  let theme = ['light', 'dark', 'system'].includes(requested) ? requested
    : ['light', 'dark', 'system'].includes(saved) ? saved : 'system';

  const icons = {
    system: 'M3 4h18v13H3z M8 21h8 M12 17v4',
    light: 'M16 12a4 4 0 1 1-8 0 4 4 0 0 1 8 0 M12 2v2 M12 20v2 M2 12h2 M20 12h2 M5 5l1.5 1.5 M17.5 17.5L19 19 M5 19l1.5-1.5 M17.5 6.5L19 5',
    dark: 'M21 13a9 9 0 0 1-10-10A9 9 0 1 0 21 13Z'
  };

  function applyTheme() {
    if (theme === 'system') delete document.documentElement.dataset.theme;
    else document.documentElement.dataset.theme = theme;
    document.querySelector('.rodin-theme-icon path')?.setAttribute('d', icons[theme]);
    document.querySelector('[aria-label="Color theme"]')?.setAttribute('title', 'Color theme: ' + theme);
    for (const option of document.querySelectorAll('.rodin-theme-menu [role="menuitemradio"]')) {
      option.setAttribute('aria-checked', String(option.dataset.theme === theme));
    }
    const light = theme === 'light' || (theme === 'system' && system.matches);
    document.querySelector('meta[name="theme-color"]')?.setAttribute('content', light ? '#f7f8fa' : '#20252d');
  }
  applyTheme();
  system.addEventListener('change', applyTheme);

  document.addEventListener('DOMContentLoaded', () => {
    const brand = document.querySelector('#m-navbar-brand');
    if (!brand) return;
    const item = document.createElement('div');
    item.className = 'rodin-theme-control';
    const trigger = document.createElement('button');
    trigger.type = 'button';
    trigger.id = 'rodin-theme-trigger';
    trigger.className = 'rodin-theme-trigger';
    trigger.setAttribute('aria-label', 'Color theme');
    trigger.setAttribute('aria-haspopup', 'menu');
    trigger.setAttribute('aria-expanded', 'false');
    trigger.setAttribute('aria-controls', 'rodin-theme-menu');
    const icon = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    icon.setAttribute('viewBox', '0 0 24 24');
    icon.setAttribute('aria-hidden', 'true');
    icon.classList.add('rodin-theme-icon');
    icon.append(document.createElementNS('http://www.w3.org/2000/svg', 'path'));
    trigger.append(icon);
    const menu = document.createElement('div');
    menu.id = 'rodin-theme-menu';
    menu.className = 'rodin-theme-menu';
    menu.setAttribute('role', 'menu');
    menu.setAttribute('aria-labelledby', trigger.id);
    menu.hidden = true;
    const options = [];

    function closeMenu(focus = false) {
      menu.hidden = true;
      trigger.setAttribute('aria-expanded', 'false');
      if (focus) trigger.focus();
    }

    function openMenu(index = options.findIndex(option => option.dataset.theme === theme)) {
      menu.hidden = false;
      trigger.setAttribute('aria-expanded', 'true');
      options[index].focus();
    }

    for (const [value, text] of [['system', 'System'], ['light', 'Light'], ['dark', 'Dark']]) {
      const option = document.createElement('button');
      option.type = 'button';
      option.setAttribute('role', 'menuitemradio');
      option.tabIndex = -1;
      option.dataset.theme = value;
      option.textContent = text;
      option.addEventListener('click', () => {
        theme = value;
        try { localStorage.setItem(key, theme); } catch {}
        const url = new URL(location.href);
        url.searchParams.delete('theme');
        history.replaceState(history.state, '', url);
        applyTheme();
        closeMenu(true);
      });
      options.push(option);
      menu.append(option);
    }
    trigger.addEventListener('click', () => menu.hidden ? openMenu() : closeMenu());
    trigger.addEventListener('keydown', event => {
      if (event.key === 'ArrowDown' || event.key === 'ArrowUp') {
        event.preventDefault();
        openMenu(event.key === 'ArrowDown' ? 0 : options.length - 1);
      }
    });
    menu.addEventListener('keydown', event => {
      const index = options.indexOf(document.activeElement);
      if (event.key === 'Escape') {
        event.preventDefault();
        closeMenu(true);
      } else if (event.key === 'Tab') {
        closeMenu(true);
      } else if (['ArrowDown', 'ArrowUp', 'Home', 'End'].includes(event.key)) {
        event.preventDefault();
        const next = event.key === 'Home' ? 0 : event.key === 'End' ? options.length - 1
          : (index + (event.key === 'ArrowDown' ? 1 : options.length - 1)) % options.length;
        options[next].focus();
      }
    });
    document.addEventListener('pointerdown', event => {
      if (!item.contains(event.target)) closeMenu();
    });
    item.addEventListener('focusout', event => {
      if (!item.contains(event.relatedTarget)) closeMenu();
    });
    item.append(trigger, menu);
    brand.closest('nav').append(item);
    applyTheme();
  });
})();
