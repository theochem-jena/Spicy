# Contributing

## Source formatting
### Documentation and Commenting
Every top level function and declaration needs to be documented. Use [Haddock](https://www.haskell.org/haddock/) to write visible documentation for every function, where assumptions and behaviour of the function are clearly described. In case there could be any ambiguities what the arguments of a function mean, also document them using Haddock. Top level functions documenting should use Haddock stly block comments with the comment opening and closing braces in extra lines.

Local declarations in `let` or `where` should be documented if it is not absolutely clear, what their meaning is. Use normal (non-Haddock) comments for this.

If you use any abbreviations, that are not immediately obvious, put a comment there, so that the function is easily understandable.

__*Example:*__
```haskell
{-|
Substitute all indices in the bonding style 'IntMap' by replacements given in a 'IntMap' 'Int'. If
keys of the old bonding map do not appear in the replacement map, a 'Left' error 'String' will be
given. Therefore the replacement map must be complete and no holes are allowed.
-}
reIndexBonds ::
     IntMap Int                    -- ^ Substitution patterns with original as key and new as Value.
  -> IntMap IntSet                 -- ^ Original '_molecule_Bonds' styled 'IntMap'.
  -> Either String (IntMap IntSet) -- ^ New bonds after replacements. Results in a 'Left' 'String'
                                   --   as error, if not all old indices can be mapped to new
                                   --   indices.
reIndexBonds repMap bondsMap
  | IS.null lostKeysIM && IS.null lostKeysIS =
      Right $ (IM.map (updateIntSet repMap)) . (updateKeys repMap) $ bondsMap
  | otherwise     =
      Left "reIndexBonds: The replacement IntMap has holes and cannot map all old indices to new \
      \ indices."
  where
    -- Old keys (origin of bonds Atoms)
    oldKeys    = IM.keysSet bondsMap
    -- Old values (target of bonds Atoms)
    oldSets    = IS.unions $ bondsMap
    -- Keys to be replaced
    repKeys    = IM.keysSet repMap
    -- Keys that are in the old map, but not in the new map.
    lostKeysIM = oldKeys \\ repKeys
    lostKeysIS = oldSets \\ repKeys
  ```
  
### Function Declaration and Indentation
Line width maximum is 100 characters, both for comments and for code. Generally functions should not be longer than 30 lines (excluding comments and possibly multiline records). If a function cannot be expressed in 30 lines, split it into atomic functions. Short functions are easier to understand, debug and reuse. Expections might be parsers and writers for text processing or wrapper interactions in a process.
  
Indentation width is exactly 2 spaces. `->` in `case` statements and `=` in assignments of the same level (e.g. in the same `where` or `let` blocks should be on the same column. The code in a `let ... in  ...` construct should begin in the same column. Therefore one space after `let` but two after `in`.

In case of multiline functions, do not start the function definition in the same line as the `=`.
  
__*Example:*__
```haskell
{-|
Checks if a replacement 'IntMap' is complete for replacing both 'IM.Key's and 'IntSet' values of
an 'IntMap' 'IntSet' and gives 'True' if no holes (no mapping from old to new value) are found.
-}
isRepMapCompleteIM ::
     IntMap Int    -- ^ Replacement 'IntMap', mapping from old 'IM.Key's to new 'IM.Key's.
  -> IntMap IntSet -- ^ Data structure to compare with.
  -> Bool          -- ^ Result. 'True' if replacement 'IntMap' is complete, 'False' otherwise.
isRepMapCompleteIM repMap bondsMap =
  let -- Old keys (origin of bonds Atoms)
      oldKeys    = IM.keysSet bondsMap
      -- Old values (target of bonds Atoms)
      oldSets    = IS.unions $ bondsMap
      -- Keys to be replaced
      repKeys    = IM.keysSet repMap
      -- Keys that are in the old map, but not in the new map.
      lostKeysIM = oldKeys \\ repKeys
      lostKeysIS = oldSets \\ repKeys
  in  if IS.null lostKeysIM && IS.null lostKeysIS
        then True
        else False
```

## Testing and Benchmarking
For all major functions (all that are exported in a module) there must be a unit or golden test. *ALL* tests muss pass. Tests are in written in the [Tasty](http://hackage.haskell.org/package/tasty) framework, both for unit tests and golden tests. 
Tests should be performed by a simple `stack test`.

Spicy shall be able to handle large molecules up to protein size. Use [Criterion](http://hackage.haskell.org/package/criterion) benchmarks for functions, that need to handle large data structures and test for different input sizes, to check if they scale and parallelise well. Benchmarks can be executed by `stack bench`, but this can cause problems with parallel functions. Better build the executable and run it standalone as `stack exec benchmarks`.
