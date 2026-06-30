// Runtime config, loaded before the app. Fill these in AFTER deploying Umami
// (see docs/DEPLOY.md → Analytics). Leave empty to run with analytics OFF.
window.DISC_ANALYTICS = {
  src: "https://umami-production-2ffd.up.railway.app/script.js",
  websiteId: "b3640a0b-b040-42da-b6bb-256c7b3ffa21",
  host: "",   // not needed: the script reports to its own origin (the Umami domain)
};
