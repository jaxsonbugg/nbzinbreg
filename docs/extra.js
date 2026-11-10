// Website-only behaviour for the report page (the PDF is not affected).
//
// The tables in the report's appendix are folded by default: each one is
// wrapped in a closed <details> element whose summary is the table's caption,
// so readers see what each table is and open the ones they want. A link to a
// folded table (for example "Table 6" in the text) opens it before scrolling
// to it. The figures in the appendix stay open.
(function () {
  "use strict";

  function appendixTables() {
    var start = document.getElementById("sec-appendix");
    if (!start) return [];
    var floats = document.querySelectorAll('div.quarto-float[id^="tbl-"]');
    var out = [];
    for (var i = 0; i < floats.length; i++) {
      // keep the tables that come after the appendix heading in the page
      if (start.compareDocumentPosition(floats[i]) & Node.DOCUMENT_POSITION_FOLLOWING) {
        out.push(floats[i]);
      }
    }
    return out;
  }

  function captionText(el) {
    var cap = el.querySelector("figcaption, caption");
    var txt = cap ? cap.textContent : el.id;
    return txt.replace(/\s+/g, " ").trim();
  }

  function fold(el) {
    if (el.closest("details.nbz-fold")) return;
    var details = document.createElement("details");
    details.className = "nbz-fold";
    var summary = document.createElement("summary");
    summary.textContent = captionText(el) + " (click to show)";
    var body = document.createElement("div");
    body.className = "nbz-fold-body";
    el.parentNode.insertBefore(details, el);
    details.appendChild(summary);
    details.appendChild(body);
    body.appendChild(el);
    details.addEventListener("toggle", function () {
      summary.textContent = captionText(el) +
        (details.open ? " (click to hide)" : " (click to show)");
    });
  }

  // Open the folded table that contains the element with this id, if any.
  function reveal(id) {
    if (!id) return;
    var target = document.getElementById(decodeURIComponent(id));
    if (!target) return;
    var details = target.closest("details.nbz-fold");
    if (details && !details.open) details.open = true;
    target.scrollIntoView();
  }

  document.addEventListener("DOMContentLoaded", function () {
    appendixTables().forEach(fold);

    document.addEventListener("click", function (event) {
      var link = event.target.closest('a[href^="#tbl-"]');
      if (!link) return;
      var id = link.getAttribute("href").slice(1);
      var target = document.getElementById(id);
      if (target && target.closest("details.nbz-fold")) {
        event.preventDefault();
        history.pushState(null, "", "#" + id);
        reveal(id);
      }
    });

    window.addEventListener("hashchange", function () {
      reveal(window.location.hash.slice(1));
    });
    if (window.location.hash) reveal(window.location.hash.slice(1));
  });
})();
