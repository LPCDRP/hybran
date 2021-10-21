
# Special Considerations

## Large Number of Reference Annotations

If you are using a large number of references, you may get an error during the annotation transfer step like this:

```
ERROR: mummer and/or mgaps returned non-zero
```

If so, you should recompile mummer using `make CPPFLAGS="-O3 -DSIXTYFOURBITS"` as described at <https://sourceforge.net/p/mummer/mailman/message/34069398/>.
