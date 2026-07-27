DotKill
=======

Spatio-temporal dotcrawl and rainbow remover for VapourSynth

Functions
=========

DotKillS
--------

dotkill.DotKillS(clip clip[, int iterations=1])

A purely spatial dotcrawl remover that can be safely used on most material after field matching.

Takes constant format YUV420 or GRAY input with any integer sample depth from 8 to 16 bits and must be at least 6x6 pixels. Only the luma plane is processed and its output is clamped to the 16-235 range scaled to the sample depth.

Note that the output of a deeper clip is not a bit exact match of the 8 bit output for the same source. Corrections are applied at full sample precision instead of being rounded to 8 bit steps, so results drift apart as iterations is raised.

iterations: The number of times to apply the internal filter. Usally a number between 1 and 4 will have the best results and using too high values may cause artifacting. Values outside the 1-10 range are clamped.

DotKillZ
--------

dotkill.DotKillZ(clip clip[, int order=0, int offset=0])

A pseudo-spatial dotcrawl and rainbow remover. It only works on NTSC content with rainbows added after 3:2 pulldown. This is true most of the time for anime.

Note that due to its nature only every other final frame will have dotcrawl and rainbows removed. Typically never artifacts if all requirements are met.

Takes constant format YUV420 input with any integer sample depth from 8 to 16 bits.

order: Field order. Usually 0, note that 1 hasn't been tested due to a lack of test material.

offset: The cycle offset for the pulldown pattern. A number between 0 and 4, other values are an error. Can only be determined by trial and error.

DotKillT
--------

dotkill.DotKillT(clip clip[, int order=0, int offset=0, int dupthresh=64, int tratio=3, bint show=False])

A full spatio-temporal dotcrawl and rainbow remover. It only works on NTSC content with rainbows added after 3:2 pulldown. This is true most of the time for anime.

May produce extreme artifacting if dupthresh is set too high.

Takes constant format YUV420 input with any integer sample depth from 8 to 16 bits.

order: Field order. Usually 0, note that 1 hasn't been tested due to a lack of test material.

offset: The cycle offset for the pulldown pattern. A number between 0 and 4, other values are an error. Can only be determined by trial and error.

dupthresh: The threshold for determining if a block has changed between fields. Depending on the source material 32-128 are usually reasonable values. A value of 0 makes the function identical to DotKillZ. Always given in 8 bit units and scaled internally for deeper input, so a given value means the same thing at any supported sample depth.

tratio: If more than 1/tratio blocks have changed between fields then temporal filtering won't be considered in that direction. Higher values can make high motion sections less likely to artifact.

show: Shows which blocks have been determined to contain no change between fields and therefore will be blended to reduce artifacts. White square means that it will blend with the next frame and black square the previous. The squares are drawn at full scale rather than in the TV range. Note that cycle offset 4 filters both fields separately and draws the squares once per field, so the two sets of squares partially overwrite each other there.

Output frame properties
-----------------------

DotKillTOffset: The cycle position the frame was processed as, meaning (n + offset) % 5. Useful for verifying that offset has been guessed correctly.

Usage
=====

clip = core.dotkill.DotKillT(clip, offset=1, dupthresh=64)

clip = core.vivtc.VFM(clip)

clip = core.dotkill.DotKillS(clip, iterations=4)

clip = core.vivtc.VDecimate(clip)

