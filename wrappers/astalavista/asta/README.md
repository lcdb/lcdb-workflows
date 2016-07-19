# Wrapper for astalavista's `asta` tool

[AStalavista home](http://sammeth.net/confluence/display/ASTA/Home)

## Input
- `gtf`: GTF file, conforming to the AStalavista requirements. Must be uncompressed.


## Output

A reformatted GTF-like file with splice events.

Note that AStalavista hard-codes the output file to
`{input.gtf}_astalavista.gtf.gz`; the wrapper will move this hard-coded file to
the output indicated in the rule.

## Threads

Does not run in parallel

## Params
* `extra`: extra arguments passed verbatim to `asta`

## Example

```python
rule asta:
    input:
    gtf="path/to.gtf"
    output: "asta.gtf"
    params: extra=" -e [ASE,ASI]"
    wrapper:
        "/path/to/wrapper/location"
```
