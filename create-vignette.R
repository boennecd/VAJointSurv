library(rmarkdown)
setwd("vignettes")

render("VAJointSurv.Rmd", output_format = html_document(
  css = "VAJointSurv.css"))
render("VAJointSurv.Rmd", output_format = github_document(
  pandoc_args = "--webtex=https://render.githubusercontent.com/render/math?math="),
  output_file = "README.md")
file.copy("README.md", file.path("..", "README.md"), overwrite = TRUE)
unlink("README.md")
unlink("README.html")

unlink(file.path("..", "man", "figures"), recursive = TRUE)
dir.create(file.path("..", "man", "figures"))
fs <- list.files(file.path("man", "figures"))
for(f in fs)
  file.copy(file.path("man", "figures", f),
            file.path("..", "man", "figures", f), overwrite = TRUE)
