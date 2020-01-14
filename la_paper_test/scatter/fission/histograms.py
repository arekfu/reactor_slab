import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

XSs = {}
labels = ["DT", "MDT", "NWDT", "MNWDT", "CT", "MCT", "ICT"]

XSs["LI"] = []
XSs["LI"].append([5.687972e-01, 1.075544e+00, 3.342886e-01, 1.933680e-01, 8.042081e-01, 3.677820e+00, 1.087019e+00])
XSs["LI"].append([3.471030e-01, 8.632787e-01, 7.424985e-01, 9.676993e+00, 5.099425e-01, 2.216417e+00, 8.702700e-01])
XSs["LI"].append([9.069280e-01, 1.720841e+00, 2.978698e-01, 2.721813e-01, 6.279362e-01, 5.531620e-01, 1.737728e+00])

XSs["LD"] = []
XSs["LD"].append([8.149060e-01, 1.168369e+00, 1.617051e-01, 1.499649e-01, 2.250754e+00, 2.067454e+01, 1.170074e+00])
XSs["LD"].append([6.055655e-01, 1.002581e+00,1.844817e+00,4.363391e+00,1.219219e+00,6.342114e+00,1.003287e+00])
XSs["LD"].append([1.302003e+00,1.861099e+00,2.253619e-01,2.270270e-01,5.057634e-01,4.354844e-01,1.872630e+00])

XSs["EI"] = []
XSs["EI"].append([1.524192e-01,5.719868e-01,1.337528e+00,2.822216e-01,1.941186e-01,8.107828e-01,5.715516e-01])
XSs["EI"].append([9.484022e-02,4.444092e-01,1.928970e-01,7.804061e-01,1.251515e-01,6.486152e-01,4.446166e-01])
XSs["EI"].append([2.053560e-01,7.672789e-01,8.755477e-02,1.820837e-01,1.967894e-01,4.876745e-01,7.684357e-01])

XSs["ED"] = []
XSs["ED"].append([2.673992e-01,8.630623e-01,1.106175e-01,2.635608e-01,2.992375e-01,9.304053e-01,8.674301e-01])
XSs["ED"].append([8.990382e-01,1.738056e+00,5.739496e-01,1.216393e+00,6.548927e-01,1.392745e+00,1.746948e+00])
XSs["ED"].append([1.813577e+01,5.883079e+01,1.952231e+00,1.617453e+01,1.175481e+01,3.224204e+01,5.784776e+01])

XSs["SG"] = []
XSs["SG"].append([4.710912e-02,1.432405e-01,1.782232e-02,5.113038e-02,4.967437e-02,1.460750e-01,1.392317e-01])
XSs["SG"].append([6.848160e-02,1.970181e-01,6.121516e-02,1.808593e-01,6.674694e-02,1.891538e-01,1.924251e-01])
XSs["SG"].append([9.681388e+01,2.839148e+02,3.804556e+00,1.817162e+01,3.290774e+01,8.993709e+01,2.907293e+02])

XSs["BG"] = []
XSs["BG"].append([9.402770e+00,1.154674e+01,2.331385e-01,2.162312e-01,2.559039e+00,1.982370e+00,1.139038e+01])
XSs["BG"].append([2.223414e+00,4.137136e+00,3.289681e+00,1.928424e+00,1.133826e+01,3.338007e+00,4.088441e+00])
XSs["BG"].append([5.009031e+00,6.110353e+00,9.353882e-01,9.240246e-01,1.431564e+00,1.300503e+00,6.112586e+00])

x = np.arange(len(labels))
width = 0.8
nbars = 3

for xsname in XSs.keys():
    fig, ax = plt.subplots()
    A = ax.bar(x-width/nbars,XSs[xsname][0],width/nbars,label="Collisions (Real) Score")
    B = ax.bar(x, XSs[xsname][1], width/nbars, label="Collisions (All) Score")
    C = ax.bar(x + width/nbars,XSs[xsname][2],width/nbars,label="Leakage Score")
    
    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                    xy=(rect.get_x()+rect.get_width() / 2, height),
                    xytext=(0,3), textcoords="offset points", ha='center',
                    va='bottom')


    #autolabel(A)
    #autolabel(B)
    #autolabel(C)
    ax.set_title(xsname + " FOM Results")

    ax.set_ylabel("Figure of Merit")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    fig.tight_layout()
    #plt.savefig(xsname+"_bar_.pgf")
    plt.show()
