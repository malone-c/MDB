demand = readxl::read_xlsx('data/MDBWaterMarketDataset_Demand_v1.0.0.xlsx', sheet = 2)
supply = readxl::read_xlsx('data/MDBWaterMarketDataset_Supply_v1.1.0.xlsx', sheet = 2)



sup = supply[supply$P != 0, ]

sup.1 = data.frame(year = sup$year, region = sup$region, 
									 price = sup$P, rainfall = sup$R)


## Get the elements of demand s.t. the year-region pair
## is included in sup.1

dem = NULL
for (i in 1:dim(sup)[1]) {
	dem.new = demand[demand$Year == sup[i,]$year & demand$Region == sup[i,]$region, ]
	dem = rbind(dem, dem.new)
}
dem

index = which(!(dem$Volume_applied %in% c("NA", "0")))

dem.na.rm = dem[index, ]
dem[dem$Volume_applied != 0,]

### There are 1616 obs in this set

df = data.frame(year = NULL,
								region = NULL,
								vol = NULL,
								price = NULL,
								crop = NULL,
								rainfall = NULL)
for (i in 1:dim(sup)[1]) {
	for (j in 1:dim(dem)[1]) {
		if (sup$year[i] == dem$Year[j]
				& sup$region[i] == dem$Region[j]
				& !(dem$Volume_applied[j] %in% c("NA", "0"))
				) {
			new.obs = data.frame(year = sup$year[i], 
													 region = sup$region[i], 
													 vol = as.numeric(dem$Volume_applied[j]),
													 price = as.numeric(sup$P[i]),
													 crop = dem$Industry[j],
													 rainfall = sup$R[i]
													 )
			df = rbind(df, new.obs)
		}
	}
}

df.cmplt = data.frame(year = NULL,
												region = NULL,
												vol = NULL,
												price = NULL,
												crop = NULL,
												rainfall = NULL)
for (i in unique(df$region)) {
	for (j in unique(df$crop[df$region == i])) {
		truth.condition = length(unique(df$year[df$region == i & df$crop == j])) == 13
		if (truth.condition) {
			subset.cmplt = df[df$region == i & df$crop == j, ]
			df.cmplt = rbind(df.cmplt, subset.cmplt)
		}
	}
}

# It is very important to understand what I have just done.
# I have gone through each crop-region pair, and I have 
# removed any pair from the dataset which did not have a full set 
# of 13 years. The reason I have done this is because missing values
# for years will cause problems with the Cochrane-Orcutt correction.

df$vol = as.numeric(df$vol)

dim(df.cmplt)
# There are 1222 observations left.

# Figure out how many observations we will get rid of in the 
# Cochrane-Orcutt correction.
count = 0
for (i in unique(df$region)) {
	for (j in unique(df$crop[df$region == i])) {
		count = count + 1
	}
}
count
# count is 152. So we will remove 152 observations.
# We will be left with 1222 - 152 = 1070 observations.

head(df)

# 159 unique price values.
length(unique(df$price))


hist(df.cmplt$vol, xlim = c(0, 400000))
df.cmplt


pairs(df.cmplt)

lm.simp = lm(vol ~ price, data = df.cmplt)
plot(lm.simp)

lm.all = lm(log(vol) ~ ., data = df.cmplt)
plot(lm.all)

plot(lm(log(vol) ~ price, data = df.cmplt), which = 2)

plot(log(vol) ~ price, data = df.cmplt)
lines(lowess(log(df.cmplt$vol) ~ df.cmplt$price), col = 'red')



summary(lm.all)

car::avPlots(lm.all, terms = "price")
anova(lm.all)

# Investigate potential serial correlation of errors
boxplot(lm.all$residuals ~ df.cmplt$year)

lines(df.cmplt$year, lowess(lm.all$residuals ~ df.cmplt$year)$y, lwd = 2, col = 'red')
abline(0, 0, col = 'red', lty = 2)

# There are some weird outlying points
plot(cooks.distance(lm.all), type = 'h')
identify(cooks.distance(lm.all))

df.cmplt[110:130, ]
# The weird points are the first 3 records for hay pastures
# in the lower MD (years 06, 07, and 08). Note that the
# volume amounts are absolutely tiny I reckon these are data
# entry errors.
summary(df.cmplt$vol)
min(df.cmplt$vol)

sum(df.cmplt$vol > 100) / length(df.cmplt$vol)
# 98.7% of observations have volumes over 100.

df.olrmvd = df.cmplt[!(df.cmplt$region == "NSW Lower Darling" & df.cmplt$crop == "Pastures - Hay"), ]
dim(df.olrmvd)
dim(df.cmplt)

lm.all.olrmvd = lm(log(vol) ~ ., data = df.olrmvd)
plot(lm.all.olrmvd)

# Residuals are left skewed

car::avPlots(lm.all.olrmvd)
	

qq = function(model) {
	qqnorm(model$residuals)
	qqline(model$residuals)
}

lm.qq = lm(I(vol^(1/3)) ~ ., data = df.olrmvd)
qq(lm.qq)

qq(lm(log(vol + 17000) ~ ., data = df.olrmvd))
	
plot(lm(log(vol) ~ ., data = df.olrmvd))
qqnorm(c(lm.qq$residuals))
qqline(c(lm.qq$residuals))
# Really nice. Cube root is a better transform than log, for
# normality.


plot(lm.qq, which = 1)
# Nonconstant variance definitely violated seriously

plot(lm.all.olrmvd, which = 1)
# So this one has a much nicer residual plot, but some left
# skew in the residuals. It is a better model than lm.qq.

anova(lm.qq)
anova(lm.all.olrmvd)

car::avPlots(lm.all.olrmvd, terms = "price")




lm.new = lm(log(vol) ~ . - year, data = df.olrmvd)
plot(lm.new)
car::avPlots(lm.new, terms = ~ price)
summary(lm.new)
summary(lm.all.olrmvd)
anova(lm.new)
#### TODO
# Remove the crop type from the model.


length(unique(df.olrmvd$price))
# 150 unique prices


### Look at an autocorrelation plot for the means

means = NULL
for (i in 1:13) {
	year = unique(df.olrmvd$year)[i]
	mean.new = mean(df.olrmvd$vol[df.olrmvd$year == year])
	means = cbind(means, mean.new)
}
means.now = means[2:13]
means.lag = means[1:12]
cor(means.now, means.lag)

# Yep


length(unique(df.olrmvd$crop))








#### TEST FOR AC

df.sort = df.olrmvd[order(df.olrmvd$region, df.olrmvd$crop, df.olrmvd$year), ]
head(df.sort)

rownames(df.sort) = seq(1, dim(df.sort)[1], 1)

n = dim(df.sort)[1]

lag <- df.lag <- NULL
for (i in seq(1, n, 13)) {
	lagged.logvol = log(df.sort$vol[i:(i+11)])
	t.rows = cbind(df.sort[(i+1):(i+12), ], lagged.logvol)
	lag = rbind(lag, lagged.logvol)
	df.lag = rbind(df.lag, t.rows)
}



lm.lag = lm(log(vol) ~ ., data = df.lag)
summary(lm.lag)
anova(lm.lag)


library(nlme)
gls = gls(log(vol) ~ ., data = df.lag)
summary(gls)
anova(gls)
plot(gls)

boxplot(log(df.lag$vol) ~ df.lag$year)
abline(mean(log(df.lag$vol)), 0, col = 'red')

-0.000980 * exp(-0.000980 *mean(df.olrmvd$price))

library(panelAR)

panel = NULL
for (i in 1:dim(df.olrmvd)[1]) {
	panel = rbind(panel, paste(df.olrmvd$region[i], df.olrmvd$crop[i]))
}
df.olrmvd = cbind(df.olrmvd, panel)

AR = panelAR(log(vol) ~ ., df.olrmvd, 
				panelVar = "panel", timeVar = "year",
				autoCorr = c("ar1"),
				panelCorrMethod = "none")

summary(AR)
plot(AR)

require(car)
durbinWatsonTest(AR$residuals,max.lag = 1)

# Not as bad as the others, close to 2 which is good.


### Time to fit an actual model now.








