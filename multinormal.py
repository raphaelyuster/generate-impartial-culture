# We consider F={1234, 4231, 3412, 3214, 2413, 1432} which we wish to rule out as generating the 4-impartial
# culture. We compute p(F,C_4) and compare it with p_{C_4} where the latter is the probability of a generated
# tournament of the 4-impartial culture to have a Hamilton cycle (equivalently, being the unique 4-vertex
# tournament with a Hamilton cycle).

from scipy.stats import multivariate_normal

orthant_box = [10 for _ in range(4)]
origin = [0, 0, 0, 0]
maxpts = 2000000000
a = 1 / 3.0

# Covariance matrix corresponding to the directed cycle (1,2,3,4,1) in the impartial cuture
cov_12341_impartial = [
    [1, -a, 0, -a],
    [-a, 1, -a, 0],
    [0, -a, 1, -a],
    [-a, 0, -a, 1]]

# compute the orthant probability
cdf_value = multivariate_normal.cdf(orthant_box, mean=origin, cov=cov_12341_impartial, lower_limit=origin,
                                    allow_singular=False, maxpts=maxpts)
print("Probability of generating the cycle  (1,2,3,4,1) =", cdf_value)
print("p_{C_4} =", cdf_value * 6)  # multiply by 6 since aut(C_4)=4

# Covariance matrix M1 corresponding to the directed cycles (1,2,3,4,1), (1,4,3,2,1) in the F-culture
cov_m1 = [
    [1, -a, a, -a],
    [-a, 1, -a, a],
    [a, -a, 1, -a],
    [-a, a, -a, 1]]
cdf_m1 = multivariate_normal.cdf(orthant_box, mean=origin, cov=cov_m1, lower_limit=origin, allow_singular=False,
                                 maxpts=maxpts)
print("Probability of generating the cycles (1,2,3,4,1) or (1,4,3,2,1) in the F-culture =", cdf_m1)

# Covariance matrix M2 corresponding to the directed cycles (1,2,4,3,1), (1,3,4,2,1) in the F-culture
cov_m2 = [
    [1, -a, -a, -a],
    [-a, 1, -a, -a],
    [-a, -a, 1, -a],
    [-a, -a, -a, 1]]
cdf_m2 = multivariate_normal.cdf(orthant_box, mean=origin, cov=cov_m2, lower_limit=origin, allow_singular=True,
                                 maxpts=maxpts)
print("Probability of generating the cycles (1,2,4,3,1) or (1,3,4,2,1) in the F-culture =", cdf_m2)

# Covariance matrix M3 corresponding to the directed cycles (1,3,2,4,1), (4,2,3,1,4) in the F-culture
cov_m3 = [
    [1, -a, a, -a],
    [-a, 1, -a, -a],
    [a, -a, 1, -a],
    [-a, -a, -a, 1]]
cdf_m3 = multivariate_normal.cdf(orthant_box, mean=origin, cov=cov_m3, lower_limit=origin, allow_singular=False,
                                 maxpts=maxpts)
print("Probability of generating the cycles (1,3,2,4,1) or (4,2,3,1,4) in the F-culture =", cdf_m3)
print("p(F,C_4) =", cdf_m1 * 2 + cdf_m2 * 2 + cdf_m3 * 2)
