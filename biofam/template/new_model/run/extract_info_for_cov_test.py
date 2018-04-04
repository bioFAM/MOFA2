import h5py

def print_attrs(name, obj):
    print name
    for key, val in obj.attrs.iteritems():
        print "    %s: %s" % (key, val)

transpose=1

if transpose:
    SigmaNode = "SigmaAlphaW"

    f = h5py.File("spatial/simul_spatial.h5", 'r')
    # f.visititems(print_attrs)
    views = f["parameters"][SigmaNode].keys()
    print("Simulated model")
    for view in views:
        print("view : " + view)
        #print("X positions cells : ", f["parameters"][SigmaNode][view]["X"][:, :])
        print("true variance hyperparameters for each factor : ", f["parameters"][SigmaNode][view]["ix"][:, ])

    f2 = h5py.File("spatial/test.hdf5", 'r')
    # f2.visititems(print_attrs)
    views = f["parameters"][SigmaNode].keys()
    print("")
    print("Fitted model")
    for view in views:
        print("view : " + view)
        print("recovered variance hyperparameters for each factor : ", f2["parameters"][SigmaNode][view]["ix"][:, ])
        print("signifiance elbo(sig=sig_recov)-elbo(sig=0) for each factor : ",
              f2["parameters"][SigmaNode][view]["sig"][:, ])

else:
    SigmaNode = "SigmaZ"
    f = h5py.File("spatial/simul_spatial.h5",'r')
    #f.visititems(print_attrs)
    print("Simulated model")
    #print("X positions cells : ",f["parameters"][SigmaNode]["X"][:,:])
    print("true variance hyperparameters for each factor : ", f["parameters"][SigmaNode]["ix"][:,])

    f2 = h5py.File("spatial/test.hdf5",'r')
    print("")
    print("Fitted model")
    print("recovered variance hyperparameters for each factor : ", f2["parameters"][SigmaNode]["ix"][:,])
    print("signifiance elbo(sig=sig_recov)-elbo(sig=0) for each factor : ", f2["parameters"][SigmaNode]["sig"][:,])
