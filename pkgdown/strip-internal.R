# Post-build fixups for the pkgdown site. Run after pkgdown::build_site(); the
# CI workflow calls this before deploying to gh-pages. Two jobs:
#
# 1. Remove CLAUDE.md. pkgdown renders every .md file in the package root
#    (pkgdown:::package_mds() has no exclusion hook), but CLAUDE.md is
#    contributor guidance for coding agents, not site content.
# 2. Add a self-referencing <link rel="canonical"> to every page. pkgdown does
#    not emit canonical tags, and GitHub Pages serves each page under more than
#    one path (e.g. `/` and `/index.html`); canonical consolidates the ranking
#    signals onto one URL per page.

site <- "site"
base_url <- "https://pedroliman.github.io/neuralsbi/"

# --- 1. Strip CLAUDE.md from HTML, sitemap, and search index -----------------

unlink(file.path(site, c("CLAUDE.html", "CLAUDE.md")))

sitemap <- file.path(site, "sitemap.xml")
x <- readLines(sitemap, warn = FALSE)
writeLines(gsub("<url><loc>[^<]*/CLAUDE\\.html</loc></url>", "", x), sitemap)

search_json <- file.path(site, "search.json")
entries <- jsonlite::read_json(search_json)
entries <- Filter(
  function(e) is.null(e$path) || !grepl("/CLAUDE\\.html", e$path),
  entries
)
jsonlite::write_json(entries, search_json, auto_unbox = TRUE, null = "null")

# --- 2. Inject canonical URLs ------------------------------------------------

html_files <- list.files(site, pattern = "\\.html$", recursive = TRUE,
                         full.names = TRUE)

for (f in html_files) {
  rel <- sub(paste0("^", site, "/"), "", f)
  # Root index.html canonicalises to the bare site URL, everything else to its
  # own path. Matches the loc entries pkgdown writes into sitemap.xml.
  href <- if (identical(rel, "index.html")) base_url else paste0(base_url, rel)
  tag <- sprintf('<link rel="canonical" href="%s">', href)

  html <- readLines(f, warn = FALSE)
  if (any(grepl("rel=\"canonical\"", html, fixed = TRUE))) next
  # Insert immediately before the first </head>.
  hpos <- grep("</head>", html, fixed = TRUE)[1]
  if (is.na(hpos)) next
  html <- append(html, tag, after = hpos - 1L)
  writeLines(html, f)
}
