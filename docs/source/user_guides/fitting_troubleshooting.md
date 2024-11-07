(Note: This content is taken directly from an old Slack post. It's possible context is missing, the contents are out of date, or the user who posted it was inexperienced with the code.)

# Fitting troubleshooting

If you spend a while trying to solve a fitting issue, then learn a simple solution, please add it here to save others time. (Also a good reference for yourself later.)

## Fit seems to be jumping up and down when it should be more smooth

Example:  
[https://imgur.com/a/R2WCR8J](https://imgur.com/a/R2WCR8J)  
**Solution:** Probably need to reduce the grid ratio of the model (`gridrat`). Basically, if the grid ratio is too high, then an image will sometimes be included and sometimes not, resulting in the jumping up and down.

## Favorite fitting does not seem to be the best fitting at all

Example:  
[https://imgur.com/a/8AxA3tR](https://imgur.com/a/8AxA3tR)  
**Solution:** Investigate if some dataset is interfering in the final result. It is possible that the data from a telescope is not as good and it is driving the fitting code to a wrong model. Compare the best model with and without that telescope. Set lower cuts for bad seeing data points from this telescope.
