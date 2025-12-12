# Web Total Station Simulator (WebTS)

A web-based application that simulates a **Leica FlexLine TS03 Total Station** for educational and field use. It is designed for manual **Stadia Tacheometry** (using a theodolite) but also supports manual Slope Distance input.

## Features

*   **Station Setup:** Known Point orientation and **Resection** (Free Station) using Least Squares Adjustment.
*   **Survey/Measure:** Calculate coordinates (E, N, H) from Stadia readings (Top/Mid/Bot wires) or Slope Distance.
*   **Stakeout:** Calculate Turn angles and Move distances to find target points.
*   **Area:** Calculate 2D area of a polygon.
*   **Data:** Save points to local storage and export as CSV.
*   **Mobile Friendly:** Works as a Progressive Web App (PWA) on phones and offline.

## How to use on your Phone (GitHub Pages)

You can host this for free on GitHub so you can access it via a URL on your phone.

1.  **Log in to GitHub.**
2.  **Create a new Repository:**
    *   Name it something like `ts03-sim`.
    *   Check "Initialize with a README" (optional).
3.  **Upload Files:**
    *   Go to your new repository page.
    *   Click **Add file** -> **Upload files**.
    *   Drag and drop ALL the files from this folder (`index.html`, `style.css`, `app.js`, `manifest.json`, `sw.js`, `icon-*.svg`).
    *   Click **Commit changes**.
4.  **Enable GitHub Pages:**
    *   Go to **Settings** (tab at the top of your repo).
    *   Scroll down to **Pages** (on the left sidebar).
    *   Under **Build and deployment** -> **Branch**, select `main` (or `master`) and folder `/ (root)`.
    *   Click **Save**.
5.  **Get your URL:**
    *   Wait about 1-2 minutes.
    *   Refresh the Pages settings page. You will see a link like `https://your-username.github.io/ts03-sim/`.
    *   Open this link on your phone!

## How to Install (Add to Home Screen)

Once you open the link on your phone:
*   **Android (Chrome):** Tap the menu (3 dots) -> **Add to Home screen** (or "Install App").
*   **iOS (Safari):** Tap the Share button -> **Add to Home Screen**.

Now it works like a native app, even without internet!
