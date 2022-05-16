# Algorithms

## Equalization of Brightness

The Equalization of Brightness algorithm normalizes the amplitude of an image without change of phases and therefore of structure. 

In general the algorithm can be applied to any image `image` which resembles a 2D array through `Eob(image)`. In this case all default values are created and can be reused. For more complex and sophisticated situations wiht mnore control, it is recommended to use the `EobAlgorithm` command, where the filterbanks and the regularitzation parameter can be stated explicitely.
To help with the creation of own steerable Riesz filterbanks, the function `EobInitialization` can be used.

```@docs
Eob
EobInitialization
EobAlgorithm
```

## multiresolution decomposition

This algorithm will follow in a future release

## multireolution orientaiton estimation

This algorithm will follow in a future release