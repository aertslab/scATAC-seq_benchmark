def calc_kde(xy):
    from scipy.stats import gaussian_kde
    return gaussian_kde(xy)(xy)