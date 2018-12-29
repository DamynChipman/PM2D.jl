"""

"""
function VortexDiffusion(panels::Array{Panel2D},
                         
                         M::Int64)

    # Extract geometry from panels
    NPAN = length(panels)
    n_hats = zeros(NPAN,2)
    t_hats = zeros(NPAN,2)
    X = zeros(NPAN,2)
    for n=1:NPAN
        n_hats[n,:] = panels[n].n_hat
        t_hats[n,:] = panels[n].t_hat
        X[n,1] = panels[n].r1[1]
        X[n,2] = panels[n].r1[2]
    end

    # Guassian Spreading and Normalization Constant
    sigma = 0.2
    a = 1.0

    function phi_ij(i,j)
        return a*exp(-norm(X(j) - X(i))/(2*sigma^2))
    end

    function grad_phi_ij(i,j)

    end

    function grad2_phi_ij(i,j)

    end





end
