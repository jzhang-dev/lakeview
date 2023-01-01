Specifying regions of interest
==============================

region: reference_name compulsory even if only one reference sequence exists for explicity

start and end coordinates are optional

samtools-compatible region string, or a two-element tuple containing reference name and a coordinate interval (start, end). If the coordinate interval is `None`, it will be assumed to be from the first base to the last base of the reference sequence