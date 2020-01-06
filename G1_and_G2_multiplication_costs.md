# Lookup table for slope an intercept

Final formula is `gas cost = intercept + slope * num_group_limbs` where `num_group_limbs` is number of 64 bits words required to represent a group order. `intercept` and `slope` are in units of gas.

## G1 (in base field)

| Number of limbs| Slope | Intercept |
|---|---|---|
|4| 855 | 4080 |
|5| 1830 | 915 |
|6| 2115 | 1275 |
|7| | | 
|8| | |
|9| | |
|10| | |
|11| | |
|12| | |
|13| | |
|14| | |
|15| | |
|16| | |

## G2 in quadratic extension

| Number of limbs| Slope | Intercept |
|---|---|---|
|4| 3345 | 12825 |
|5| 6675 | 3210 |
|6| 8325 | 1980 |
|7| | | 
|8| | |
|9| | |
|10| | |
|11| | |
|12| | |
|13| | |
|14| | |
|15| | |
|16| | |

## G2 in cubic extension

| Number of limbs| Slope | Intercept |
|---|---|---|
|4| 6825 | 17175 |
|5| 11715 | 3210 |
|6| 15240 | 2595 |
|7| | | 
|8| | |
|9| | |
|10| | |
|11| | |
|12| | |
|13| | |
|14| | |
|15| | |
|16| | |