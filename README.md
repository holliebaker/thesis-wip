# Hollie's thesis

Based on [rstudio/bookdown-demo](https://github.com/rstudio/bookdown-demo).

Serve a preview:
```R
# first time you may need to install:
install.packages("bookdown")

# then run a live preview
bookdown::serve_book(dir = ".", output_dir = "_book", preview = TRUE, in_session = TRUE, quiet = FALSE)

# or publish
bookdown::render_book()
```

Then go to `http://127.0.0.1:4321` to see a live preview ([docs](https://bookdown.org/yihui/bookdown/serve-the-book.html)).

