Changelog
**************

Version 1.1.0
==============

Reorganisation of .vtu saving by splitting into saving mesh fields and data fields. 
Mesh fields are created only for one frame and copied to the remaining frames before the time domain integration loop. 
The data fields are appended to existing frames inside the time domain integration loop.
The Matlab base64 encode was replaced by a more efficient Apache version.

Version 1.0.0
==============

Base version.
