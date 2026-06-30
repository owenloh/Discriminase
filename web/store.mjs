// Browser persistence:
//  - IndexedDB holds built indexes (the big typed arrays) so a panel is built once.
//  - localStorage holds small settings: saved panels + the last-used parameters,
//    so an accidental reload doesn't lose your typed-up set of bacteria.

const DB_NAME = "discriminase";
const STORE = "indexes";

function openDB() {
  return new Promise((resolve, reject) => {
    const req = indexedDB.open(DB_NAME, 1);
    req.onupgradeneeded = () => req.result.createObjectStore(STORE);
    req.onsuccess = () => resolve(req.result);
    req.onerror = () => reject(req.error);
  });
}

export async function saveIndex(key, payload) {
  const db = await openDB();
  return new Promise((resolve, reject) => {
    const tx = db.transaction(STORE, "readwrite");
    tx.objectStore(STORE).put(payload, key);
    tx.oncomplete = () => resolve();
    tx.onerror = () => reject(tx.error);
  });
}

export async function loadIndex(key) {
  const db = await openDB();
  return new Promise((resolve, reject) => {
    const tx = db.transaction(STORE, "readonly");
    const rq = tx.objectStore(STORE).get(key);
    rq.onsuccess = () => resolve(rq.result || null);
    rq.onerror = () => reject(rq.error);
  });
}

export async function listIndexKeys() {
  const db = await openDB();
  return new Promise((resolve, reject) => {
    const tx = db.transaction(STORE, "readonly");
    const rq = tx.objectStore(STORE).getAllKeys();
    rq.onsuccess = () => resolve(rq.result);
    rq.onerror = () => reject(rq.error);
  });
}

// --- small settings (localStorage) ----------------------------------------
export function saveSetting(key, value) {
  try { localStorage.setItem("disc:" + key, JSON.stringify(value)); } catch (e) { /* ignore quota */ }
}
export function loadSetting(key, fallback) {
  try { const s = localStorage.getItem("disc:" + key); return s ? JSON.parse(s) : fallback; }
  catch (e) { return fallback; }
}
