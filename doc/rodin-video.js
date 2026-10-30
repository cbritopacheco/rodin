/* Player controls adapted from https://cbritopacheco.github.io/assets/js/video-motion.js. */
(() => {
  const reducedMotion = window.matchMedia("(prefers-reduced-motion: reduce)");
  const touchInput = window.matchMedia("(hover: none), (pointer: coarse)");
  const icons = {
    play: '<path d="m9 5 11 7-11 7Z" fill="currentColor" stroke="none"/>',
    pause: '<path d="M8 5v14M16 5v14" stroke-width="4"/>',
    muted: '<path d="M11 5 6 9H3v6h3l5 4ZM16 9l6 6m0-6-6 6"/>',
    volume: '<path d="M11 5 6 9H3v6h3l5 4ZM16 8a6 6 0 0 1 0 8m3-11a10 10 0 0 1 0 14"/>',
    expand: '<path d="M8 3H3v5m13-5h5v5M3 16v5h5m13-5v5h-5"/>',
    collapse: '<path d="M3 8h5V3m8 0v5h5M8 21v-5H3m13 5v-5h5"/>'
  };
  function icon(name) {
    return `<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" focusable="false">${icons[name]}</svg>`;
  }
  function time(seconds) {
    if (!Number.isFinite(seconds)) return "0:00";
    return `${Math.floor(seconds / 60)}:${String(Math.floor(seconds % 60)).padStart(2, "0")}`;
  }
  const players = [...document.querySelectorAll("main article video")].map((video) => {
    const frame = document.createElement("div");
    frame.className = "rodin-video-player";
    frame.classList.toggle("uses-touch-controls", touchInput.matches);
    touchInput.addEventListener("change", () => {
      frame.classList.toggle("uses-touch-controls", touchInput.matches);
      frame.classList.remove("controls-visible");
    });
    frame.style.maxWidth = video.style.maxWidth;
    frame.style.margin = video.style.margin;
    video.style.margin = "0";
    video.before(frame);
    frame.append(video);
    const controls = document.createElement("div");
    controls.className = "rodin-video-controls";
    controls.setAttribute("role", "group");
    controls.setAttribute("aria-label", `${video.getAttribute("aria-label") || "Video"} controls`);
    function button(label, symbol) {
      const el = document.createElement("button");
      el.type = "button";
      el.setAttribute("aria-label", label);
      el.title = label;
      el.innerHTML = icon(symbol);
      controls.append(el);
      return el;
    }
    const play = button("Pause video", "pause");
    const clock = document.createElement("span");
    clock.className = "rodin-video-time";
    clock.setAttribute("aria-hidden", "true");
    clock.textContent = "0:00";
    controls.append(clock);
    const seek = document.createElement("input");
    seek.type = "range";
    seek.min = "0";
    seek.max = "1000";
    seek.value = "0";
    seek.setAttribute("aria-label", "Video position");
    controls.append(seek);
    const mute = button("Unmute video", "muted");
    const fullscreen = button("Enter fullscreen", "expand");
    if (!frame.requestFullscreen && !video.webkitEnterFullscreen) fullscreen.hidden = true;
    frame.append(controls);
    // Keep player navigation independent of the documentation search shortcuts.
    frame.addEventListener("keydown", (event) => event.stopPropagation());
    video.controls = false;
    video.tabIndex = 0;
    function toggle() {
      if (video.paused) video.play().catch(() => {});
      else video.pause();
    }
    function update() {
      frame.classList.toggle("is-paused", video.paused);
      play.innerHTML = icon(video.paused ? "play" : "pause");
      play.setAttribute("aria-label", video.paused ? "Play video" : "Pause video");
      play.title = play.getAttribute("aria-label");
      mute.innerHTML = icon(video.muted ? "muted" : "volume");
      mute.setAttribute("aria-label", video.muted ? "Unmute video" : "Mute video");
      mute.title = mute.getAttribute("aria-label");
    }
    function updateTime() {
      const duration = video.duration;
      seek.disabled = !Number.isFinite(duration) || duration <= 0;
      if (!seek.disabled) seek.value = String(Math.round(1000 * video.currentTime / duration));
      seek.setAttribute("aria-valuetext", `${time(video.currentTime)} of ${time(duration)}`);
      clock.textContent = time(video.currentTime);
    }
    play.addEventListener("click", toggle);
    frame.addEventListener("pointerdown", (event) => {
      if (event.pointerType === "touch") frame.classList.add("uses-touch-controls");
      else frame.classList.remove("controls-visible");
    });
    video.addEventListener("click", () => {
      if (frame.classList.contains("uses-touch-controls")) frame.classList.toggle("controls-visible");
      else toggle();
    });
    frame.addEventListener("focusin", (event) => {
      if (event.target.matches(":focus-visible")) frame.classList.add("controls-visible");
    });
    frame.addEventListener("focusout", (event) => {
      if (!frame.contains(event.relatedTarget)) frame.classList.remove("controls-visible");
    });
    document.addEventListener("pointerdown", (event) => {
      if (!frame.contains(event.target)) frame.classList.remove("controls-visible");
    });
    video.addEventListener("keydown", (event) => {
      if (event.key === " " || event.key === "Enter") { event.preventDefault(); toggle(); }
    });
    mute.addEventListener("click", () => { video.muted = !video.muted; });
    seek.addEventListener("input", () => {
      if (Number.isFinite(video.duration)) video.currentTime = Number(seek.value) / 1000 * video.duration;
    });
    fullscreen.addEventListener("click", () => {
      if (document.fullscreenElement === frame) document.exitFullscreen().catch(() => {});
      else if (frame.requestFullscreen) frame.requestFullscreen().catch(() => {});
      else if (video.webkitEnterFullscreen) video.webkitEnterFullscreen();
    });
    document.addEventListener("fullscreenchange", () => {
      const active = document.fullscreenElement === frame;
      fullscreen.innerHTML = icon(active ? "collapse" : "expand");
      fullscreen.setAttribute("aria-label", active ? "Exit fullscreen" : "Enter fullscreen");
      fullscreen.title = fullscreen.getAttribute("aria-label");
    });
    ["play", "pause", "ended", "volumechange"].forEach((event) => video.addEventListener(event, update));
    ["timeupdate", "loadedmetadata", "durationchange"].forEach((event) => video.addEventListener(event, updateTime));
    update();
    updateTime();
    return video;
  });
  // Preserve each tutorial's playback policy; pause when reduced motion is requested.
  function motion() {
    if (reducedMotion.matches) players.forEach((video) => video.pause());
  }
  motion();
  reducedMotion.addEventListener("change", motion);
})();
