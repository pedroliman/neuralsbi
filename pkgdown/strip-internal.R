# Remove CLAUDE.md from the rendered site. pkgdown renders every .md file in
# the package root (pkgdown:::package_mds() has no exclusion hook), but
# CLAUDE.md is contributor guidance for coding agents, not site content.
# Run after pkgdown::build_site(); the CI workflow calls this before deploying.
site <- "site"

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
