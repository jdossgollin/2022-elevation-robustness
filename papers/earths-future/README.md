## Earth's Future submission

To generate the document tracking changes from `submission_01` to `submission_02`:

```shell
latexdiff papers/earths-future/submission_01.tex papers/earths-future/submission_02.tex > papers/earths-future/trackedchanges_01_02.tex
```

To fix problems with compiling the resulting document, see [this StackOverflow answer](https://tex.stackexchange.com/questions/574280/latexdiff-with-cite-commands-gives-output-with-apparently-mismatched-braces/579662#579662).
It involves editing the `latexdiff` file which requires `sudo` permissions, but it's worked here!
