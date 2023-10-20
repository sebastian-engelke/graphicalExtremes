# Three arguments in chunk hooks

A chunk hook has three arguments: `before`, `options` and `envir`. We show how they work through some simple examples.

## The `before` argument

It is a logical argument: `before == TRUE` executes code before a chunk.


```r
library(knitr)
now <- Sys.time()
knit_hooks$set(foo1 = function(before) {
  if (before) {
    '_I appear before a chunk!_\n\n'
  } else {
    '\n\n_I am after a chunk..._'
  }
})
now <- Sys.time()
knit_hooks$set(timeit = function(before, options, envir) {
    if (before) {
        now <<- Sys.time()
        deparse(options)
    } else {
        paste(sprintf("\nChunk rendering time: %s seconds.\n", round(Sys.time() - now, digits = 3))) 
    }
})

knitr::opts_chunk$set(
    collapse = TRUE,
    comment = '#>',
    error = TRUE,
    timeit = TRUE
)
```


structure(list(eval = TRUE, echo = TRUE, results = "markup",     tidy = FALSE, tidy.opts = NULL, collapse = TRUE, prompt = FALSE,     comment = "#>", highlight = TRUE, size = "normalsize", background = "#F7F7F7",     strip.white = FALSE, cache = 0, cache.path = "cache/", cache.vars = NULL,     cache.lazy = TRUE, dependson = NULL, autodep = FALSE, cache.rebuild = FALSE,     fig.keep = "high", fig.show = "asis", fig.align = "default",     fig.path = "figure/", dev = "png", dev.args = NULL, dpi = 72,     fig.ext = NULL, fig.width = 7, fig.height = 7, fig.env = "figure",     fig.cap = NULL, fig.scap = NULL, fig.lp = "fig:", fig.subcap = NULL,     fig.pos = "", out.width = NULL, out.height = NULL, out.extra = NULL,     fig.retina = 1, external = TRUE, sanitize = FALSE, interval = 1,     aniopts = "controls,loop", warning = TRUE, error = TRUE,     message = TRUE, render = NULL, ref.label = NULL, child = NULL,     engine = "R", split = FALSE, include = TRUE, purl = TRUE,     timeit = TRUE, label = "unnamed-chunk-2", code = c("library(knitr)",     "now <- Sys.time()", "print('set hook?')"), out.width.px = 504,     out.height.px = 504, params.src = ""), class = "knitr_strict_list")

```r
library(knitr)
now <- Sys.time()
print('set hook?')
#> [1] "set hook?"
```


Chunk rendering time: 0.005 seconds.


structure(list(eval = TRUE, echo = TRUE, results = "markup",     tidy = FALSE, tidy.opts = NULL, collapse = TRUE, prompt = FALSE,     comment = "#>", highlight = TRUE, size = "normalsize", background = "#F7F7F7",     strip.white = FALSE, cache = 0, cache.path = "cache/", cache.vars = NULL,     cache.lazy = TRUE, dependson = NULL, autodep = FALSE, cache.rebuild = FALSE,     fig.keep = "high", fig.show = "asis", fig.align = "default",     fig.path = "figure/", dev = "png", dev.args = NULL, dpi = 72,     fig.ext = NULL, fig.width = 7, fig.height = 7, fig.env = "figure",     fig.cap = NULL, fig.scap = NULL, fig.lp = "fig:", fig.subcap = NULL,     fig.pos = "", out.width = NULL, out.height = NULL, out.extra = NULL,     fig.retina = 1, external = TRUE, sanitize = FALSE, interval = 1,     aniopts = "controls,loop", warning = TRUE, error = TRUE,     message = TRUE, render = NULL, ref.label = NULL, child = NULL,     engine = "R", split = FALSE, include = TRUE, purl = TRUE,     timeit = FALSE, label = "unnamed-chunk-3", code = c("x <- 100",     "print(x)"), out.width.px = 504, out.height.px = 504, params.src = "timeit=FALSE"), class = "knitr_strict_list")

```r
x <- 100
print(x)
#> [1] 100
```


Chunk rendering time: 0.004 seconds.


Test the `foo1` hook:

structure(list(eval = TRUE, echo = TRUE, results = "markup",     tidy = FALSE, tidy.opts = NULL, collapse = TRUE, prompt = FALSE,     comment = "#>", highlight = TRUE, size = "normalsize", background = "#F7F7F7",     strip.white = FALSE, cache = 0, cache.path = "cache/", cache.vars = NULL,     cache.lazy = TRUE, dependson = NULL, autodep = FALSE, cache.rebuild = FALSE,     fig.keep = "high", fig.show = "asis", fig.align = "default",     fig.path = "figure/", dev = "png", dev.args = NULL, dpi = 72,     fig.ext = NULL, fig.width = 7, fig.height = 7, fig.env = "figure",     fig.cap = NULL, fig.scap = NULL, fig.lp = "fig:", fig.subcap = NULL,     fig.pos = "", out.width = NULL, out.height = NULL, out.extra = NULL,     fig.retina = 1, external = TRUE, sanitize = FALSE, interval = 1,     aniopts = "controls,loop", warning = TRUE, error = TRUE,     message = TRUE, render = NULL, ref.label = NULL, child = NULL,     engine = "R", split = FALSE, include = TRUE, purl = TRUE,     timeit = TRUE, foo1 = "whatever", label = "unnamed-chunk-4",     code = "1+1", out.width.px = 504, out.height.px = 504, params.src = "foo1='whatever'"), class = "knitr_strict_list")_I appear before a chunk!_

```r
1+1
#> [1] 2
```



_I am after a chunk..._
Chunk rendering time: 0.004 seconds.
