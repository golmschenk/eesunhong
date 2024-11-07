(Note: This content is taken directly from an old Slack post. It's possible context is missing, the contents are out of date, or the user who posted it was inexperienced with the code.)

# Data source considerations and what to do with them

Note, these are not hard rules. They are just common steps you’ll often need to perform, but sometimes they might not be appropriate.

## KMT data

KMT preparation:
- Probably remove any data points with very large flux error (perhaps anything greater than `2000`).
- Probably remove any data points with terrible seeing (`FWHM = -1`).
- May also need to remove data points with bad seeing (`FWHM > 10` is getting pretty bad).

## MOA data

MOA data from Ian’s reduction:
The MOA data provided by Ian already has data points with bad seeing and large flux error removed. The only necessary preparation step may be to subtract 2450000.0 from the columns for the days.
(TODO)

## microFUN data

Amateur astronomer data that can often help fill the gaps, provide missing critical information, or otherwise augment the primary datasets.
(TODO)
