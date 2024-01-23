## Operator Design

### The `HString`
`HString`- `Λ×16` `i32`(alias for `Int32`) Array. Row indices are operator orders and column indices are for the internal storage of one operator.
```
struct HString
    buf::Vector{i32}
end
```
functions like getindex and set index need to be implemented.

### The `Operator`
One `Operator` is a column slice of length 16 of the HString. The 16 `i32` are for:

```
 0 + (1:5) | self info, 28 bits for j and 4 bits for ψ
 5 + (1:5) | prev legs, 28 bits for p and 4 bits for r
10 + (1:5) | next legs, 28 bits for p and 4 bits for r
16         |      flag, 28 bits for p and 4 bits for type
```

### Update scheme

