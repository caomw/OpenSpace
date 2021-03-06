--[[  OpenSpace keybinding script ]]--

-- Load the common helper functions
dofile(openspace.absPath('${SCRIPTS}/common.lua'))

openspace.clearKeys()
helper.setCommonKeys()
helper.setDeltaTimeKeys({
--  1           2           3           4           5           6           7           8           9           0
--------------------------------------------------------------------------------------------------------------------------
--  1s          2s          5s          10s         30s         1m          2m          5m          10m         30m
    1,          2,          5,          10,         30,         60,         120,        300,        600,        1800,

--  1h          2h          3h          6h          12h         1d          2d          4d          1w          2w
    3600,       7200,       10800,      21600,      43200,      86400,      172800,     345600,     604800,     1209600,

--  1mo         2mo         3mo         6mo         1yr         2y          5y          10y         20y         50y
    2592000,    5184000,    7776000,    15552000,   31536000,   63072000,   157680000,  315360000,  630720000,  1576800000
})
--  OBS: One month (1mo) is approximated by 30 days.
