import gwsurrogate

def install_surrogates():
    """Download the NRSur7dq4 and NRHybSur3dq8 surrogate data via gwsurrogate."""
    print('Installing NRSur7dq4...')
    gwsurrogate.catalog.pull('NRSur7dq4')

    print('NRSur7dq4 installed. Installing NRHybSur3dq8...')
    gwsurrogate.catalog.pull('NRHybSur3dq8')

    print('NRHybSur3dq8 installed.')
