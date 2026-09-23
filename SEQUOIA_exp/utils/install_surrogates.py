import gwsurrogate

def install_surrogates():

    print('Installing NRSur7dq4...')
    gwsurrogate.catalog.pull('NRSur7dq4')

    print('NRSur7dq4 installed. Installing NRHybSur3dq8...')
    gwsurrogate.catalog.pull('NRHybSur3dq8')

    print('NRHybSur3dq8 installed.')