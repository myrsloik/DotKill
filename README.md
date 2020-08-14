DotKill
=======

Spatio-temporal dotcrawl and rainbow remover for VapourSynth

Functions
=========

DotKillS
--------

dotkill.DotKillS(clip clip[, int iterations=1, bint usematch=False])

A purely spatial dotcrawl remover that can be safely used on most material after field matching.

iterations: The number of times to apply the internal filter. Usally a number between 1 and 4 will have the best results and using too high values may cause artifacting.

usematch: If true then matching hints from VFM are used when processing. This may or may not have a positive effect.

DotKillZ
--------

dotkill.DotKillZ(clip clip[, int order=0, int offset=0])

A pseudo-spatial dotcrawl and rainbow remover. It only workes on NTSC content with rainbows added after 3:2 pulldown. This is true most of the time for anime.

Note that due to its nature only every other final frame will have dotcrawl and rainbows removed. Typically never artifacts if all requirements are met.

order: Field order. Usually 0, note that 1 hasn't been tested due to a lack of test material.

offset: The cycle offset for the pulldown pattern. A number between 0 and 4. Can only be determined by trial and error.

DotKillT
--------

dotkill.DotKillT(clip clip[, int order=0, int offset=0, int dupthresh=64, int tratio=3, bint show=False])

A full spatioi-temporal dotcrawl and rainbow remover. It only workes on NTSC content with rainbows added after 3:2 pulldown. This is true most of the time for anime.

May produce extreme artifacting if dupthresh is set too high.

order: Field order. Usually 0, note that 1 hasn't been tested due to a lack of test material.

offset: The cycle offset for the pulldown pattern. A number between 0 and 4. Can only be determined by trial and error.

dupthresh: The threshold for determining if a block has changed between fields. Depending on the source material 32-128 are usually reasonable values. A value of 0 makes the function identical to DotKillZ.

tratio: If more than 1/tratio blocks have changed between fields then temporal filtering won't be considered in that direction. Higher values can make high motion sections less likely to artifact.

show: Shows which blocks have been determined to contain no change between frames and therefore will be blended to reduce artifacts. White square means that it will blend with the next frame and black square the previous. 

Usage
=====

clip = core.dotkill.DotKillT(clip, offset=1, dupthresh=64)
clip = core.vivtc.VFM(clip)
clip = core.dotkill.DotKillS(clip, iterations=4)
clip = core.vivtc.VDecimate(clip)

