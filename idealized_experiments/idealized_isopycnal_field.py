import numpy as np


class densityField:
    def __init__(self, g=10, Nsquared=1e-5, alphax=1.1e-3, alphay=1e-3, kappax=2e-5/np.pi, kappay=2e-5/np.pi, rho0=1025):
        self.g = g                # gravitational acceleration [ms^-2]
        self.Nsquared = Nsquared  # Square of buoyancy frequency [s^-2]
        self.alphax = alphax      # zonal amplitude
        self.alphay = alphay      # meridional amplitude
        self.kappax = kappax      # zonal wavenumber [m^-1]
        self.kappay = kappay      # meridional wavenumber [m^-1]
        self.rho0 = rho0           # background density [kg m^-3]

        
    def create_interpolated_density_grid(self, nx=101, ny=201, nz=301, Lx=1e6, Ly=2e6, H=3000, rhotype='averaged'):
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.Lx = Lx
        self.Ly = Ly
        self.H = H
        self.X = np.linspace(0, Lx, nx)
        self.Y = np.linspace(0, Ly, ny)
        self.Z = np.linspace(0, -H, nz)

        self.dX = np.abs(self.X[1] - self.X[0])
        self.dY = np.abs(self.Y[1] - self.Y[0])
        self.dZ = np.abs(self.Z[1] - self.Z[0])

        self.XXX, self.YYY, self.ZZZ = np.meshgrid(self.X, self.Y, self.Z, indexing='ij')

        if rhotype == 'interpolated':
            rho_interpolated = self.rho0 * (
                1
                - self.Nsquared * self.ZZZ / self.g
                + self.alphax * np.sin(self.kappax * self.XXX)
                + self.alphay * np.sin(self.kappay * self.YYY)
            )
            self.rho_interpolated
    
        if rhotype == 'averaged':
            # Create staggered versions of original grid.
            self.create_staggered_grid()
            integrated_rho = self.rho0 * (self.XXX_stag * self.YYY_stag * self.ZZZ_stag 
                             - self.Nsquared * self.ZZZ_stag **2 * self.XXX_stag * self.YYY_stag / (2 * self.g) 
                             - self.alphax * np.cos(self.kappax * self.XXX_stag) * self.YYY_stag * self.ZZZ_stag / self.kappax 
                             - self.alphay * np.cos(self.kappay * self.YYY_stag) * self.XXX_stag * self.ZZZ_stag / self.kappay)
            self.rho_averaged = np.diff(np.diff(np.diff(integrated_rho, axis=0), axis=1), axis=2)/(self.dX * self.dY * self.dZ)
            
            
    def create_staggered_grid(self):
        if not hasattr(self, "XXX"):
            raise RuntimeError("Cannot create a stagered grid before initializing an interpolated density grid (use `create_interpolated_density_grid()` first).")
        self.X_stag = np.linspace(-self.dX, self.Lx, self.nx+1) + self.dX/2
        self.Y_stag = np.linspace(-self.dY, self.Ly, self.ny+1) + self.dY/2
        # Due to the definition of Z in Parcels!
        self.Z_stag = np.linspace(self.dZ, -self.H, self.nz+1) - self.dZ/2

        self.XXX_stag, self.YYY_stag, self.ZZZ_stag = np.meshgrid(self.X_stag, self.Y_stag, self.Z_stag, indexing='ij')
        
        
    def compute_slope(self):
        if not hasattr(self, "XXX"):
            raise RuntimeError("Cannot compute a staggered rho field before initializing an interpolated density grid (use `create_interpolated_density_grid()` first).")
        if not hasattr(self, 'X_stag'):
            self.create_staggered_grid()
        
        XXX_stagX, YYY_stagX, ZZZ_stagX = np.meshgrid(self.X_stag, self.Y, self.Z, indexing='ij')
        XXX_stagY, YYY_stagY, ZZZ_stagY = np.meshgrid(self.X, self.Y_stag, self.Z, indexing='ij')
        XXX_stagZ, YYY_stagZ, ZZZ_stagZ = np.meshgrid(self.X, self.Y, self.Z_stag, indexing='ij')
        
        rho_stagX = self.rho0 * (
                1
                - self.Nsquared * ZZZ_stagX / self.g
                + self.alphax * np.sin(self.kappax * XXX_stagX)
                + self.alphay * np.sin(self.kappay * YYY_stagX)
            )
        
        rho_stagY = self.rho0 * (
                1
                - self.Nsquared * ZZZ_stagY / self.g
                + self.alphax * np.sin(self.kappax * XXX_stagY)
                + self.alphay * np.sin(self.kappay * YYY_stagY)
            )
        
        rho_stagZ = self.rho0 * (
                1
                - self.Nsquared * ZZZ_stagZ / self.g
                + self.alphax * np.sin(self.kappax * XXX_stagZ)
                + self.alphay * np.sin(self.kappay * YYY_stagZ)
            )
        
        drho_dx = np.diff(rho_stagX, axis=0)/self.dX # np.gradient(rho, axis=1)
        drho_dy = np.diff(rho_stagY, axis=1)/self.dY # np.gradient(rho, axis=1)
        drho_dz = np.diff(rho_stagZ, axis=2)/self.dZ # np.gradient(rho, axis=0)
        
        # Due to the definition depth in Parcels! Normally these should be negative
        self.Sx = drho_dx/drho_dz
        self.Sy = drho_dy/drho_dz
        self.Sabs2 = self.Sx**2 + self.Sy**2
        
    def compute_slope_deriv(self):
        if not hasattr(self, "XXX"):
            raise RuntimeError("Cannot compute a staggered rho field before initializing an interpolated density grid (use `create_interpolated_density_grid()` first).")
        if not hasattr(self, 'X_stag'):
            self.create_staggered_grid()
        if not hasattr(self, 'Sx'):
            self.compute_slope()
        self.dSxdx = np.diff(self.Sx, axis=0)/self.dX
        self.dSxdy = np.diff(self.Sx, axis=1)/self.dY
        self.dSxdz = np.diff(self.Sx, axis=2)/self.dZ
        self.dSydx = np.diff(self.Sy, axis=0)/self.dX
        self.dSydy = np.diff(self.Sy, axis=1)/self.dY
        self.dSydz = np.diff(self.Sy, axis=2)/self.dZ
        

    def isopycnal_array(self, x, y, rho_iso=1040):
        z_iso = self.g / self.Nsquared* (1 - rho_iso / self.rho0 + self.alphax * np.sin(self.kappax * x)+ self.alphay * np.sin(self.kappay * y))
        return z_iso
    

    def isopycnal_grid(self, rho, x_iso, y_iso):
        XX, YY = np.meshgrid(x_iso, y_iso, indexing='ij')
        ZZ = (
            self.g
            / self.Nsquared
            * (
                (1 - rho / self.rho0)
                + self.alphax * np.sin(self.kappax * XX)
                + self.alphay * np.sin(self.kappay * YY)
            )
        )
        return ZZ, XX, YY