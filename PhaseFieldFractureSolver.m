% MATLAB Class File
% PhaseFieldFractureSolver - AT2 Phase Field Fracture Solver for Anti-Plane Shear
% 
% This is a MATLAB class implementing the AT2 phase field fracture model.
% For more details, see the documentation and references within.

% MIT License
%
% Copyright (c) 2026 Yongxing Shen
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

classdef PhaseFieldFractureSolver
    % PhaseFieldFractureSolver - AT2 Phase Field Fracture Solver
    %
    %   This class implements the AT2 phase field fracture model for simulating
    %   anti-plane shear (Mode III) fracture behavior in brittle materials.
    %
    %   Theoretical Background:
    %       - Degradation function: g(d) = (1-d)^2
    %       - Fracture energy: psi_frac = Gc/l0 * d^2 + Gc*l0 * |grad(d)|^2
    %       - History variable: H = max(psi_elastic) for irreversibility
    %
    %   References:
    %       沈泳星. 变分断裂相场法基础[M]. 北京: 清华大学出版社, 2026.
    %       SHEN Y X. Variational Phase Field Method for Fracture: Fundamentals[M]. Beijing: Tsinghua University Press, 2026. (in Chinese)
    %
    %   Example:
    %       mesh = PhaseFieldFractureSolver.createMesh(51, 51, 1.0, 1.0);
    %       mat = PhaseFieldFractureSolver.createMaterial(210e9, 0.3, 1e1, 0.02);
    %       [u, d] = PhaseFieldFractureSolver.solveStaggered(mesh, mat, u0, d0, ...);

    %===========================================================================
    %                        MESH CLASS
    %===========================================================================
    methods (Static)
        function mesh = createMesh(nx, ny, Lx, Ly)
            % createMesh - Create a rectangular mesh for 2D analysis
            %   mesh = createMesh(nx, ny, Lx, Ly) creates a structured mesh
            %   with nx nodes in x-direction, ny nodes in y-direction, covering
            %   a domain of size Lx x Ly.
            %
            %   Input:
            %       nx - number of nodes in x-direction (integer >= 2)
            %       ny - number of nodes in y-direction (integer >= 2)
            %       Lx - domain length in x-direction (meters)
            %       Ly - domain length in y-direction (meters)
            %
            %   Output:
            %       mesh - structure containing mesh data with fields:
            %           .nx, .ny - node counts
            %           .Lx, .Ly - domain dimensions
            %           .dx, .dy - element dimensions
            %           .x, .y - node coordinates
            %           .elements - element connectivity (nElements x 4)
            %           .nNodes, .nElements, .nDOF, .nPhase
            %           .H, .Hprev - history variables
            %
            %   Example:
            %       mesh = PhaseFieldFractureSolver.createMesh(51, 51, 1.0, 1.0);
            mesh = struct();
            mesh.nx = int32(nx);
            mesh.ny = int32(ny);
            mesh.Lx = Lx;
            mesh.Ly = Ly;
            mesh.dx = Lx / (nx - 1);
            mesh.dy = Ly / (ny - 1);

            mesh.x = zeros(nx, ny);
            mesh.y = zeros(nx, ny);
            [mesh.x, mesh.y] = meshgrid(linspace(0, Lx, nx), linspace(0, Ly, ny));
            mesh.x = mesh.x';
            mesh.y = mesh.y';

            mesh.nNodes = nx * ny;
            mesh.nElements = (nx-1) * (ny-1);
            mesh.nDOF = int32(mesh.nNodes);
            mesh.nPhase = int32(mesh.nNodes);

            mesh.elements = zeros(mesh.nElements, 4, 'int32');
            e = 1;
            for j = 1:ny-1
                for i = 1:nx-1
                    n1 = (j-1)*nx + i;
                    n2 = (j-1)*nx + i + 1;
                    n3 = j*nx + i + 1;
                    n4 = j*nx + i;
                    mesh.elements(e,:) = [n1, n2, n3, n4];
                    e = e + 1;
                end
            end

            mesh.H = zeros(nx-1, ny-1, 4);
            mesh.Hprev = zeros(nx-1, ny-1, 4);
            mesh.currentLoadStep = int32(0);
        end

        function mesh = startNewLoadStep(mesh)
            % startNewLoadStep - Advance to next load step
            %   mesh = startNewLoadStep(mesh) copies current history variable H to
            %   Hprev and increments the load step counter. This prepares the mesh
            %   for the next load increment in quasi-static analysis.
            %
            %   See also: createMesh, solveStaggered
            mesh.Hprev = mesh.H;
            mesh.currentLoadStep = mesh.currentLoadStep + 1;
        end

        function showMesh(mesh)
            % showMesh - Display mesh properties
            %   showMesh(mesh) prints a formatted summary of mesh properties
            %   including grid size, domain dimensions, element sizes, and DOF counts.
            %
            %   See also: createMesh
            fprintf('Mesh Properties:\n');
            fprintf('  Grid Size: %d x %d nodes\n', mesh.nx, mesh.ny);
            fprintf('  Domain Size: Lx = %.4f m, Ly = %.4f m\n', mesh.Lx, mesh.Ly);
            fprintf('  Element Size: dx = %.6f m, dy = %.6f m\n', mesh.dx, mesh.dy);
            fprintf('  Number of Elements: %d\n', mesh.nElements);
            fprintf('  Number of Nodes: %d\n', mesh.nNodes);
            fprintf('  Number of DOFs (displacement): %d\n', mesh.nDOF);
            fprintf('  Number of Phase Field DOFs: %d\n', mesh.nPhase);
        end
    end

    %===========================================================================
    %                        MATERIAL CLASS
    %===========================================================================
    methods (Static)
        function mat = createMaterial(E, nu, Gc, l0)
            % createMaterial - Create material properties structure
            %   mat = createMaterial(E, nu, Gc, l0) defines material properties
            %   for the AT2 phase field fracture model.
            %
            %   Input:
            %       E   - Young's modulus (Pa)
            %       nu  - Poisson's ratio (-)
            %       Gc  - Critical energy release rate (J/m^2)
            %       l0  - Length scale parameter (m)
            %
            %   Output:
            %       mat - structure with fields: E, nu, Gc, l0, G
            %
            %   Note: Shear modulus G is computed as G = E / (2*(1+nu))
            %
            %   Example:
            %       mat = PhaseFieldFractureSolver.createMaterial(210e9, 0.3, 1e1, 0.02);
            %
            %   See also: showMaterial
            mat = struct();
            mat.E = E;
            mat.nu = nu;
            mat.Gc = Gc;
            mat.l0 = l0;
            mat.G = E / (2 * (1 + nu));
        end

        function showMaterial(mat)
            % showMaterial - Display material properties
            %   showMaterial(mat) prints a formatted summary of material properties
            %   including elastic modulus, Poisson's ratio, shear modulus,
            %   fracture energy, and length scale.
            %
            %   See also: createMaterial
            fprintf('Material Properties:\n');
            fprintf('  Young''s Modulus E = %.2e Pa\n', mat.E);
            fprintf('  Poisson''s Ratio nu = %.3f\n', mat.nu);
            fprintf('  Shear Modulus G = %.2e Pa\n', mat.G);
            fprintf('  Critical Energy Release Rate Gc = %.2e J/m^2\n', mat.Gc);
            fprintf('  Length Scale Parameter l0 = %.4f m\n', mat.l0);
            fprintf('  Fracture Model: AT2\n');
        end
    end

    %===========================================================================
    %                    ELEMENT COMPUTATIONS
    %===========================================================================
    properties (Constant, Hidden)
        GAUSS_PTS = [-1/sqrt(3), 1/sqrt(3)];
        GAUSS_WTS = [1.0, 1.0];
    end

    methods (Static)
        function N = shapeFunctions(xi, eta)
            % shapeFunctions - Compute bilinear shape functions at given point
            %   N = shapeFunctions(xi, eta) returns the 4 shape function values
            %   at a point with natural coordinates (xi, eta) in a 4-node quad element.
            %
            %   Input:
            %       xi  - natural coordinate in x-direction (-1 to 1)
            %       eta - natural coordinate in y-direction (-1 to 1)
            %
            %   Output:
            %       N   - 1x4 array of shape function values
            %
            %   Shape functions:
            %       N1 = 0.25*(1-xi)*(1-eta)
            %       N2 = 0.25*(1+xi)*(1-eta)
            %       N3 = 0.25*(1+xi)*(1+eta)
            %       N4 = 0.25*(1-xi)*(1+eta)
            N = 0.25 * [ (1-xi).*(1-eta), (1+xi).*(1-eta), ...
                         (1+xi).*(1+eta), (1-xi).*(1+eta) ];
        end

        function [dN_dxi, dN_deta] = shapeFunctionDerivatives(xi, eta)
            % shapeFunctionDerivatives - Compute derivatives of shape functions
            %   [dN_dxi, dN_deta] = shapeFunctionDerivatives(xi, eta) returns the
            %   partial derivatives of shape functions with respect to natural
            %   coordinates at point (xi, eta).
            %
            %   Input:
            %       xi  - natural coordinate in x-direction (-1 to 1)
            %       eta - natural coordinate in y-direction (-1 to 1)
            %
            %   Output:
            %       dN_dxi   - derivatives w.r.t. xi (1x4)
            %       dN_deta  - derivatives w.r.t. eta (1x4)
            dN_dxi = 0.25 * [ -(1-eta), (1-eta), (1+eta), -(1+eta) ];
            dN_deta = 0.25 * [ -(1-xi), -(1+xi), (1+xi), (1-xi) ];
        end

        function [J, detJ, invJ, dN_dx, dN_dy] = computeJacobian(mesh, elem_nodes, dN_dxi, dN_deta)
            % computeJacobian - Compute Jacobian transformation from natural to physical coordinates
            %   [J, detJ, invJ, dN_dx, dN_dy] = computeJacobian(mesh, elem_nodes, dN_dxi, dN_deta)
            %   computes the Jacobian matrix and its inverse, along with the shape
            %   function derivatives in physical coordinates.
            %
            %   Input:
            %       mesh       - mesh structure
            %       elem_nodes - array of 4 node indices for the element
            %       dN_dxi    - shape function derivatives w.r.t. xi
            %       dN_deta   - shape function derivatives w.r.t. eta
            %
            %   Output:
            %       J      - 2x2 Jacobian matrix
            %       detJ   - determinant of J
            %       invJ   - inverse of J
            %       dN_dx  - shape function derivatives w.r.t. x
            %       dN_dy  - shape function derivatives w.r.t. y
            coords = [mesh.x(elem_nodes)', mesh.y(elem_nodes)'];
            J = [dN_dxi; dN_deta] * coords;
            detJ = J(1,1) * J(2,2) - J(1,2) * J(2,1);
            if abs(detJ) < 1e-15
                detJ = sign(detJ) * 1e-15;
            end
            invJ = (1/detJ) * [J(2,2), -J(1,2); -J(2,1), J(1,1)];
            dN_dx = invJ(1,1) * dN_dxi + invJ(1,2) * dN_deta;
            dN_dy = invJ(2,1) * dN_dxi + invJ(2,2) * dN_deta;
        end

        function B = computeBMatrix(dN_dx, dN_dy)
            % computeBMatrix - Compute strain-displacement matrix for anti-plane shear
            %   B = computeBMatrix(dN_dx, dN_dy) forms the strain-displacement
            %   matrix for anti-plane (Mode III) analysis.
            %
            %   For anti-plane shear with displacement w(x,y):
            %       strain = [gamma_xz, gamma_yz]^T = [dw/dx, dw/dy]^T
            %
            %   Input:
            %       dN_dx - shape function derivatives w.r.t. x (1x4)
            %       dN_dy - shape function derivatives w.r.t. y (1x4)
            %
            %   Output:
            %       B - 2x4 strain-displacement matrix
            B = [dN_dx; dN_dy];
        end

        function G = computePhaseGradientMatrix(dN_dx, dN_dy)
            % computePhaseGradientMatrix - Compute phase field gradient matrix
            %   G = computePhaseGradientMatrix(dN_dx, dN_dy) forms the gradient
            %   operator for the phase field variable d.
            %
            %   Input:
            %       dN_dx - shape function derivatives w.r.t. x
            %       dN_dy - shape function derivatives w.r.t. y
            %
            %   Output:
            %       G - 2x4 gradient matrix
            G = [dN_dx; dN_dy];
        end

        function [gp_idx, w] = getGaussIndexAndWeight(gp1, gp2)
            % getGaussIndexAndWeight - Get Gaussian quadrature index and weight
            %   [gp_idx, w] = getGaussIndexAndWeight(gp1, gp2) returns the
            %   Gaussian point index and combined weight for 2D integration.
            %
            %   Input:
            %       gp1 - first Gaussian point coordinate
            %       gp2 - second Gaussian point coordinate
            %
            %   Output:
            %       gp_idx - combined Gaussian point index (1-4)
            %       w      - combined weight
            idx1 = 1 + (gp1 >= 0);
            idx2 = 1 + (gp2 >= 0);
            w = PhaseFieldFractureSolver.GAUSS_WTS(idx1) * PhaseFieldFractureSolver.GAUSS_WTS(idx2);
            gp_idx = (idx1 - 1) * 2 + idx2;
        end

        function [Ke, Re] = computeDisplacementElementContribution(mesh, mat, e, u, d)
            % computeDisplacementElementContribution - Compute element stiffness and residual for displacement
            %   [Ke, Re] = computeDisplacementElementContribution(mesh, mat, e, u, d)
            %   computes the element stiffness matrix Ke and residual vector Re for
            %   the displacement field in the AT2 phase field model.
            %
            %   The weak form of the equilibrium equation:
            %       integral(grad(w) . sigma) dV = 0
            %   where sigma = g(d) * G * strain
            %
            %   Input:
            %       mesh - mesh structure
            %       mat  - material properties
            %       e    - element number
            %       u    - global displacement vector
            %       d    - global phase field vector
            %
            %   Output:
            %       Ke - 4x4 element stiffness matrix
            %       Re - 4x1 element residual vector
            %
            %   See also: assembleDisplacementStiffness, assembleDisplacementResidual
            nodes = mesh.elements(e, :);
            Ke = zeros(4, 4);
            Re = zeros(4, 1);
            nodeU = nodes;

            for gp_i = 1:2
                for gp_j = 1:2
                    xi = PhaseFieldFractureSolver.GAUSS_PTS(gp_i);
                    eta = PhaseFieldFractureSolver.GAUSS_PTS(gp_j);
                    N = PhaseFieldFractureSolver.shapeFunctions(xi, eta);
                    [dN_dxi, dN_deta] = PhaseFieldFractureSolver.shapeFunctionDerivatives(xi, eta);
                    [~, detJ, ~, dN_dx, dN_dy] = PhaseFieldFractureSolver.computeJacobian(mesh, nodes, dN_dxi, dN_deta);
                    B = PhaseFieldFractureSolver.computeBMatrix(dN_dx, dN_dy);

                    u_elem = u(nodeU);
                    strain = B * u_elem;
                    d_val = dot(N, d(nodes));
                    g = PhaseFieldFractureSolver.degradationFunction(d_val);
                    D = g * mat.G * eye(2);
                    sigma = D * strain;

                    [~, w] = PhaseFieldFractureSolver.getGaussIndexAndWeight(xi, eta);
                    Re = Re + B' * sigma * detJ * w;
                    Ke = Ke + B' * D * B * detJ * w;
                end
            end
        end

        function [Ke, Re] = computePhaseElementContribution(mesh, mat, e, u, d)
            % computePhaseElementContribution - Compute element stiffness and residual for phase field
            %   [Ke, Re] = computePhaseElementContribution(mesh, mat, e, u, d) computes
            %   the element stiffness matrix Ke and residual vector Re for the phase field.
            %
            %   The phase field equation comes from the variational formulation:
            %       (g'(d)*H*phi + Gc/l0*d*phi + Gc*l0*grad(d) . grad(phi)) dV = 0
            %
            %   Input:
            %       mesh - mesh structure
            %       mat  - material properties
            %       e    - element number
            %       u    - displacement vector
            %       d    - phase field vector
            %
            %   Output:
            %       Ke - 4x4 element stiffness matrix
            %       Re - 4x1 element residual vector
            %
            %   See also: assemblePhaseStiffness, assemblePhaseResidual
            nodes = mesh.elements(e, :);
            Ke = zeros(4, 4);
            Re = zeros(4, 1);
            nodeD = nodes;

            elem_i = mod(double(e)-1, double(mesh.nx)-1) + 1;
            elem_j = floor((double(e)-1) / (double(mesh.nx)-1)) + 1;

            for gp_i = 1:2
                for gp_j = 1:2
                    xi = PhaseFieldFractureSolver.GAUSS_PTS(gp_i);
                    eta = PhaseFieldFractureSolver.GAUSS_PTS(gp_j);
                    N = PhaseFieldFractureSolver.shapeFunctions(xi, eta);
                    [dN_dxi, dN_deta] = PhaseFieldFractureSolver.shapeFunctionDerivatives(xi, eta);
                    [~, detJ, ~, dN_dx, dN_dy] = PhaseFieldFractureSolver.computeJacobian(mesh, nodes, dN_dxi, dN_deta);
                    G_mat = PhaseFieldFractureSolver.computePhaseGradientMatrix(dN_dx, dN_dy);

                    d_val = dot(N, d(nodeD));
                    g_double_prime = 2.0;
                    g_prime = -2 * (1 - d_val);

                    nodeU = nodes;
                    B = PhaseFieldFractureSolver.computeBMatrix(dN_dx, dN_dy);
                    u_elem = u(nodeU);
                    strain = B * u_elem;
                    psi_t = 0.5 * mat.G * max(0, strain(1)^2 + strain(2)^2);

                    [gp_idx, w] = PhaseFieldFractureSolver.getGaussIndexAndWeight(xi, eta);
                    H_val = mesh.H(elem_i, elem_j, gp_idx);

                    Re = Re + (g_prime * N' * H_val + mat.Gc / mat.l0 * N' * d_val + mat.Gc * mat.l0 * G_mat' * (G_mat * d(nodeD))) * detJ * w;
                    Ke = Ke + (g_double_prime * H_val * (N' * N) + mat.Gc / mat.l0 * (N' * N) + mat.Gc * mat.l0 * G_mat' * G_mat) * detJ * w;
                end
            end
        end

        function g = degradationFunction(d)
            % degradationFunction - AT2 degradation function
            %   g = degradationFunction(d) computes the degradation function
            %   for the AT2 phase field model.
            %
            %   Formula: g(d) = (1-d)^2 + eps
            %   where eps = 1e-10 is a small value to ensure positive definiteness
            %
            %   Input:
            %       d - phase field value (0 = intact, 1 = fully broken)
            %
            %   Output:
            %       g - degradation function value
            %
            %   Note: The degradation function reduces stiffness as damage increases:
            %       d=0 (no damage): g = 1 (full stiffness)
            %       d=1 (fully broken): g = 0 (no stiffness)
            g = (1 - d)^2 + 1e-10;
        end

        function energy = computeElementEnergy(mesh, mat, e, u, d)
            % computeElementEnergy - Compute elastic strain energy at element Gauss points
            %   energy = computeElementEnergy(mesh, mat, e, u, d) computes the
            %   elastic strain energy density at all Gauss points of element e.
            %
            %   Energy density: psi = 0.5 * g(d) * G * |epsilon|^2
            %
            %   Input:
            %       mesh - mesh structure
            %       mat  - material properties
            %       e    - element number
            %       u    - global displacement vector
            %       d    - global phase field vector
            %
            %   Output:
            %       energy - total elastic energy at element Gauss points
            %
            %   See also: updateHistoryVariable
            nodes = mesh.elements(e, :);
            energy = 0;
            nodeD = nodes;
            nodeU = nodes;

            for gp_i = 1:2
                for gp_j = 1:2
                    xi = PhaseFieldFractureSolver.GAUSS_PTS(gp_i);
                    eta = PhaseFieldFractureSolver.GAUSS_PTS(gp_j);
                    N = PhaseFieldFractureSolver.shapeFunctions(xi, eta);
                    [dN_dxi, dN_deta] = PhaseFieldFractureSolver.shapeFunctionDerivatives(xi, eta);
                    [~, detJ, ~, dN_dx, dN_dy] = PhaseFieldFractureSolver.computeJacobian(mesh, nodes, dN_dxi, dN_deta);

                    B = PhaseFieldFractureSolver.computeBMatrix(dN_dx, dN_dy);
                    G_mat = PhaseFieldFractureSolver.computePhaseGradientMatrix(dN_dx, dN_dy);

                    u_elem = u(nodeU);
                    strain = B * u_elem;
                    d_val = dot(N, d(nodeD));

                    g = PhaseFieldFractureSolver.degradationFunction(d_val);
                    sigma = mat.G * strain;
                    psi_elastic = 0.5 * g * dot(strain, sigma);

                    d_grad = G_mat * d(nodeD);
                    psi_frac = 0.5 * mat.Gc / mat.l0 * d_val^2 + ...
                               0.5 * mat.Gc * mat.l0 * dot(d_grad, d_grad);

                    [~, w] = PhaseFieldFractureSolver.getGaussIndexAndWeight(xi, eta);
                    energy = energy + (psi_elastic + psi_frac) * detJ * w;
                end
            end
        end

        function mesh = updateHistoryVariable(mesh, mat, e, u, d)
            % updateHistoryVariable - Update history variable H at a single Gauss point
            %   mesh = updateHistoryVariable(mesh, mat, e, u, d) updates the history
            %   variable H for element e. The history variable stores the maximum elastic
            %   strain energy density ever reached, ensuring crack irreversibility.
            %
            %   Update rule: H_new = max(H_prev, psi_t)
            %   where psi_t = 0.5 * G * |epsilon|^2 is the current strain energy density
            %
            %   Input:
            %       mesh - mesh structure with H and Hprev arrays
            %       mat  - material properties
            %       e    - element number
            %       u    - global displacement vector
            %       d    - global phase field vector
            %
            %   Output:
            %       mesh - updated mesh structure
            %
            %   See also: updateHistoryVariables
            nodes = mesh.elements(e, :);
            nodeU = nodes;

            elem_i = mod(double(e)-1, double(mesh.nx)-1) + 1;
            elem_j = floor((double(e)-1) / (double(mesh.nx)-1)) + 1;

            for gp_i = 1:2
                for gp_j = 1:2
                    xi = PhaseFieldFractureSolver.GAUSS_PTS(gp_i);
                    eta = PhaseFieldFractureSolver.GAUSS_PTS(gp_j);
                    N = PhaseFieldFractureSolver.shapeFunctions(xi, eta);
                    [dN_dxi, dN_deta] = PhaseFieldFractureSolver.shapeFunctionDerivatives(xi, eta);
                    [~, detJ, ~, dN_dx, dN_dy] = PhaseFieldFractureSolver.computeJacobian(mesh, nodes, dN_dxi, dN_deta);

                    B = PhaseFieldFractureSolver.computeBMatrix(dN_dx, dN_dy);
                    u_elem = u(nodeU);
                    strain = B * u_elem;
                    psi_t = 0.5 * mat.G * max(0, strain(1)^2 + strain(2)^2);

                    [gp_idx, ~] = PhaseFieldFractureSolver.getGaussIndexAndWeight(xi, eta);
                    % 正确的历史变量更新规则：H = max(H_prev, psi_t)
                    % 使用上一步的Hprev值进行比较，确保不可逆性
                    if psi_t > mesh.Hprev(elem_i, elem_j, gp_idx)
                        mesh.H(elem_i, elem_j, gp_idx) = psi_t;
                    end
                end
            end
        end
    end

    %===========================================================================
    %                        GLOBAL ASSEMBLY
    %===========================================================================
    methods (Static)
        function K = assembleDisplacementStiffness(mesh, mat, u, d)
            % assembleDisplacementStiffness - Assemble global displacement stiffness matrix
            %   K = assembleDisplacementStiffness(mesh, mat, u, d) assembles the global
            %   stiffness matrix for the displacement field by looping over all elements.
            %
            %   Input:
            %       mesh - mesh structure
            %       mat  - material properties
            %       u    - displacement vector
            %       d    - phase field vector
            %
            %   Output:
            %       K - sparse nDOF x nDOF stiffness matrix
            %
            %   See also: assembleDisplacementResidual, computeDisplacementElementContribution
            nDOF = double(mesh.nDOF);
            K = sparse(nDOF, nDOF);

            for e = 1:mesh.nElements
                [Ke, ~] = PhaseFieldFractureSolver.computeDisplacementElementContribution(mesh, mat, e, u, d);
                nodeU = mesh.elements(e, :);
                for a = 1:4
                    for b = 1:4
                        K(nodeU(a), nodeU(b)) = K(nodeU(a), nodeU(b)) + Ke(a, b);
                    end
                end
            end
        end

        function R = assembleDisplacementResidual(mesh, mat, u, d)
            % assembleDisplacementResidual - Assemble displacement residual vector
            %   R = assembleDisplacementResidual(mesh, mat, u, d) assembles the global
            %   residual vector for the displacement equilibrium equation.
            %
            %   Input:
            %       mesh - mesh structure
            %       mat  - material properties
            %       u    - displacement vector
            %       d    - phase field vector
            %
            %   Output:
            %       R - nDOF x 1 residual vector
            %
            %   See also: assembleDisplacementStiffness
            nDOF = double(mesh.nDOF);
            R = zeros(nDOF, 1);

            for e = 1:mesh.nElements
                [~, Re] = PhaseFieldFractureSolver.computeDisplacementElementContribution(mesh, mat, e, u, d);
                nodeU = mesh.elements(e, :);
                for a = 1:4
                    R(nodeU(a)) = R(nodeU(a)) + Re(a);
                end
            end
        end

        function K = assemblePhaseStiffness(mesh, mat, u, d)
            % assemblePhaseStiffness - Assemble phase field stiffness matrix
            %   K = assemblePhaseStiffness(mesh, mat, u, d) assembles the global
            %   stiffness matrix for the phase field equation.
            %
            %   Input:
            %       mesh - mesh structure
            %       mat  - material properties
            %       u    - displacement vector
            %       d    - phase field vector
            %
            %   Output:
            %       K - sparse nPhase x nPhase stiffness matrix
            %
            %   See also: assemblePhaseResidual, computePhaseElementContribution
            nPhase = double(mesh.nPhase);
            K = sparse(nPhase, nPhase);

            for e = 1:mesh.nElements
                [Ke, ~] = PhaseFieldFractureSolver.computePhaseElementContribution(mesh, mat, e, u, d);
                nodeD = mesh.elements(e, :);
                for a = 1:4
                    for b = 1:4
                        K(nodeD(a), nodeD(b)) = K(nodeD(a), nodeD(b)) + Ke(a, b);
                    end
                end
            end
        end

        function R = assemblePhaseResidual(mesh, mat, u, d)
            % assemblePhaseResidual - Assemble phase field residual vector
            %   R = assemblePhaseResidual(mesh, mat, u, d) assembles the global
            %   residual vector for the phase field evolution equation.
            %
            %   Input:
            %       mesh - mesh structure
            %       mat  - material properties
            %       u    - displacement vector
            %       d    - phase field vector
            %
            %   Output:
            %       R - nPhase x 1 residual vector
            %
            %   See also: assemblePhaseStiffness
            nPhase = double(mesh.nPhase);
            R = zeros(nPhase, 1);

            for e = 1:mesh.nElements
                [~, Re] = PhaseFieldFractureSolver.computePhaseElementContribution(mesh, mat, e, u, d);
                nodeD = mesh.elements(e, :);
                for a = 1:4
                    R(nodeD(a)) = R(nodeD(a)) + Re(a);
                end
            end
        end

        function energy = computeTotalEnergy(mesh, mat, u, d)
            % computeTotalEnergy - Compute total elastic strain energy
            %   energy = computeTotalEnergy(mesh, mat, u, d) computes the total
            %   elastic strain energy by summing over all elements.
            %
            %   Input:
            %       mesh - mesh structure
            %       mat  - material properties
            %       u    - displacement vector
            %       d    - phase field vector
            %
            %   Output:
            %       energy - total elastic energy
            %
            %   See also: computeElementEnergy
            energy = 0;
            for e = 1:mesh.nElements
                energy = energy + PhaseFieldFractureSolver.computeElementEnergy(mesh, mat, e, u, d);
            end
        end

        function mesh = updateHistoryVariables(mesh, mat, u, d)
            % updateHistoryVariables - Update history variable for all elements
            %   mesh = updateHistoryVariables(mesh, mat, u, d) loops through all elements
            %   and updates the history variable H for each element.
            %
            %   See also: updateHistoryVariable
            for e = 1:mesh.nElements
                mesh = PhaseFieldFractureSolver.updateHistoryVariable(mesh, mat, e, u, d);
            end
        end
    end

    %===========================================================================
    %                        SOLVERS
    %===========================================================================
    methods (Static)
        function [u, history, converged, msg] = solveDisplacement(mesh, mat, d, u0, ...
                fixedDOFs, prescribedDOFs, prescribedValues, tolerance, maxIterations)
            % solveDisplacement - Newton-Raphson solver for displacement field
            %   [u, history, converged, msg] = solveDisplacement(mesh, mat, d, u0, ...)
            %   solves the displacement equilibrium equation using Newton-Raphson
            %   iteration with the phase field d held fixed.
            %
            %   The displacement equilibrium: K(u) * u = R(u,d)
            %
            %   Input:
            %       mesh              - mesh structure
            %       mat               - material properties
            %       d                 - phase field vector (fixed during displacement solve)
            %       u0                - initial displacement guess
            %       fixedDOFs         - DOFs with zero displacement (boundary)
            %       prescribedDOFs    - DOFs with prescribed displacement
            %       prescribedValues  - values of prescribed displacements
            %       tolerance         - convergence tolerance
            %       maxIterations     - maximum Newton iterations
            %
            %   Output:
            %       u         - converged displacement solution
            %       history   - convergence history structure
            %       converged - boolean indicating convergence
            %       msg       - convergence message
            %
            %   See also: solveStaggered, solvePhaseField

            u = u0;
            history = struct('iter', [], 'residual', [], 'error', []);
            converged = false;
            msg = '';

            allDOFs = 1:double(mesh.nDOF);
            freeDOFs = setdiff(allDOFs, [fixedDOFs, prescribedDOFs]);

            u_with_bc = PhaseFieldFractureSolver.applyBoundaryConditions(u, fixedDOFs, prescribedDOFs, prescribedValues);
            initialResidual = norm(PhaseFieldFractureSolver.assembleDisplacementResidual(mesh, mat, u_with_bc, d));
            initialResidual = max(initialResidual, 1e-12);
            tol = initialResidual * 1e-6;
            tol = max(tol, tolerance);

            for iter = 1:maxIterations
                u_full = PhaseFieldFractureSolver.applyBoundaryConditions(u, fixedDOFs, prescribedDOFs, prescribedValues);

                K = PhaseFieldFractureSolver.assembleDisplacementStiffness(mesh, mat, u_full, d);
                R = PhaseFieldFractureSolver.assembleDisplacementResidual(mesh, mat, u_full, d);

                K_ff = K(freeDOFs, freeDOFs);
                R_f = R(freeDOFs);

                residualNorm = norm(R_f);

                history(iter).iter = iter;
                history(iter).residual = residualNorm;

                fprintf('    Displacement Newton iter %d: ||r|| = %.6e\n', iter, residualNorm);

                if residualNorm < tol
                    converged = true;
                    msg = 'Converged';
                    fprintf('    Displacement Newton converged in %d iterations\n', iter);
                    u = PhaseFieldFractureSolver.applyBoundaryConditions(u, fixedDOFs, prescribedDOFs, prescribedValues);
                    return;
                end

                if residualNorm > 1e12
                    msg = 'Residual overflow';
                    fprintf('    CONVERGENCE FAILURE: Residual too large\n');
                    return;
                end

                try
                    delta_u_f = -K_ff \ R_f;
                catch
                    msg = 'Singular matrix';
                    fprintf('    CONVERGENCE FAILURE: Singular matrix\n');
                    return;
                end

                if any(isnan(delta_u_f)) || any(isinf(delta_u_f))
                    msg = 'NaN or Inf in increment';
                    fprintf('    CONVERGENCE FAILURE: NaN or Inf detected\n');
                    return;
                end

                u(freeDOFs) = u(freeDOFs) + delta_u_f;
            end

            msg = 'Max iterations reached';
            fprintf('    Displacement Newton: Max iterations (%d) reached\n', maxIterations);
        end

        function [d, history, converged, msg] = solvePhaseField(mesh, mat, u, d0, tolerance, maxIterations)
            % solvePhaseField - Newton-Raphson solver for phase field
            %   [d, history, converged, msg] = solvePhaseField(mesh, mat, u, d0, tolerance, maxIterations)
            %   solves the phase field evolution equation using Newton-Raphson iteration
            %   with the displacement field u held fixed.
            %
            %   The phase field equation: K_d(d) * d = R_d(u,d)
            %
            %   Input:
            %       mesh          - mesh structure
            %       mat           - material properties
            %       u             - displacement vector (fixed)
            %       d0            - initial phase field guess
            %       tolerance     - convergence tolerance
            %       maxIterations - maximum Newton iterations
            %
            %   Output:
            %       d         - converged phase field solution
            %       history   - convergence history structure
            %       converged - boolean indicating convergence
            %       msg       - convergence message
            %
            %   See also: solveStaggered, solveDisplacement
            d = d0;
            history = struct('iter', [], 'residual', []);
            converged = false;
            msg = '';

            for iter = 1:maxIterations
                K = PhaseFieldFractureSolver.assemblePhaseStiffness(mesh, mat, u, d);
                R = PhaseFieldFractureSolver.assemblePhaseResidual(mesh, mat, u, d);

                residualNorm = norm(R);

                history(iter).iter = iter;
                history(iter).residual = residualNorm;

                fprintf('    Phase field Newton iter %d: ||r|| = %.6e\n', iter, residualNorm);

                if residualNorm < tolerance
                    converged = true;
                    msg = 'Converged';
                    fprintf('    Phase field Newton converged in %d iterations\n', iter);
                    return;
                end

                if residualNorm > 1e12
                    msg = 'Residual overflow';
                    fprintf('    PHASE FIELD CONVERGENCE FAILURE: Residual too large\n');
                    return;
                end

                delta_d = -K \ R;

                if any(isnan(delta_d)) || any(isinf(delta_d))
                    msg = 'NaN or Inf in increment';
                    fprintf('    PHASE FIELD CONVERGENCE FAILURE: NaN or Inf detected\n');
                    return;
                end

                d = d + delta_d;

                fprintf('    Phase field: min=%.6f, max=%.6f, mean=%.6f\n', min(d), max(d), mean(d));
            end

            msg = 'Max iterations reached';
            fprintf('    Phase field Newton: Max iterations (%d) reached\n', maxIterations);
        end

        function [u, d, mesh, history, converged, msg] = solveStaggered(mesh, mat, u0, d0, ...
                fixedDOFs, prescribedDOFs, prescribedValues, tolerance, maxIterations)
            % solveStaggered - Staggered Newton-Raphson solver for coupled problem
            %   [u, d, mesh, history, converged, msg] = solveStaggered(mesh, mat, u0, d0, ...)
            %   solves the coupled displacement-phase field problem using a staggered
            %   (alternating) solution scheme.
            %
            %   Algorithm:
            %       1. Solve displacement with fixed d (Newton-Raphson)
            %       2. Update history variables H
            %       3. Solve phase field with fixed u (Newton-Raphson)
            %       4. Check convergence based on solution changes
            %       5. Repeat until convergence or max iterations
            %
            %   Input:
            %       mesh              - mesh structure
            %       mat               - material properties
            %       u0                - initial displacement
            %       d0                - initial phase field
            %       fixedDOFs         - DOFs with zero displacement
            %       prescribedDOFs    - DOFs with prescribed displacement
            %       prescribedValues  - values of prescribed displacements
            %       tolerance         - convergence tolerance
            %       maxIterations     - maximum staggered iterations
            %
            %   Output:
            %       u         - converged displacement solution
            %       d         - converged phase field solution
            %       mesh      - mesh with updated history variables
            %       history   - convergence history structure
            %       converged - boolean indicating convergence
            %       msg       - convergence message
            %
            %   Example:
            %       [u, d, mesh, history, converged] = PhaseFieldFractureSolver.solveStaggered(...
            %           mesh, mat, u0, d0, fixedDOFs, prescribedDOFs, values, 1e-6, 100);
            %
            %   See also: solveDisplacement, solvePhaseField
            u = u0;
            d = d0;

            history = struct('iter', [], 'uError', [], 'dError', [], ...
                'uResidual', [], 'dResidual', [], 'energy', [], ...
                'phaseConverged', [], 'dispConverged', []);
            converged = false;
            msg = '';

            u_bc = PhaseFieldFractureSolver.applyBoundaryConditions(u, fixedDOFs, prescribedDOFs, prescribedValues);
            initialUResidual = norm(PhaseFieldFractureSolver.assembleDisplacementResidual(mesh, mat, u_bc, d));
            initialUResidual = max(initialUResidual, 1e-12);
            
            % 使用相对容差和绝对容差的组合
            u_tol = max(initialUResidual * 5e-2, tolerance);
            d_tol = tolerance;
            
            % 位移和相场变化的收敛容差
            delta_u_tol = tolerance;
            delta_d_tol = tolerance;

            for iter = 1:maxIterations
                fprintf('\n=== Staggered iteration %d ===\n', iter);

                u_prev = u;
                d_prev = d;

                fprintf('Step 1: Solving displacement field...\n');
                [u, dispHistory, dispConverged, dispMsg] = ...
                    PhaseFieldFractureSolver.solveDisplacement(mesh, mat, d, u, ...
                    fixedDOFs, prescribedDOFs, prescribedValues, tolerance, 200);

                mesh = PhaseFieldFractureSolver.updateHistoryVariables(mesh, mat, u, d);

                u_bc = PhaseFieldFractureSolver.applyBoundaryConditions(u, fixedDOFs, prescribedDOFs, prescribedValues);
                currentUResidual = norm(PhaseFieldFractureSolver.assembleDisplacementResidual(mesh, mat, u_bc, d));
                phaseResidual = norm(PhaseFieldFractureSolver.assemblePhaseResidual(mesh, mat, u, d));

                if iter == 1
                    d_tol = max(phaseResidual * 1e-3, 1e-6);
                end

                err_u = norm(u - u_prev);
                    err_d = norm(d - d_prev);
                    
                    if iter > 1 && (err_u < delta_u_tol || err_d < delta_d_tol)
                        % 当位移或相场变化很小时，认为已经收敛到固定点
                        totalEnergy = PhaseFieldFractureSolver.computeTotalEnergy(mesh, mat, u, d);

                        history(iter).iter = iter;
                        history(iter).uError = err_u;
                        history(iter).dError = err_d;
                        history(iter).uResidual = currentUResidual;
                        history(iter).dResidual = phaseResidual;
                        history(iter).energy = totalEnergy;
                        history(iter).phaseConverged = true;
                        history(iter).dispConverged = dispConverged;

                        converged = true;
                        msg = sprintf('Converged in %d iterations (fixed point reached)', iter);
                        fprintf('Fixed point reached: ||Delta u|| = %.6e, ||Delta d|| = %.6e\n', err_u, err_d);
                        fprintf('\n=== STAGGERED SOLVER CONVERGED (fixed point) ===\n');
                        return;
                    end
                    
                    if phaseResidual < d_tol && dispConverged && iter > 1
                        totalEnergy = PhaseFieldFractureSolver.computeTotalEnergy(mesh, mat, u, d);

                        history(iter).iter = iter;
                        history(iter).uError = err_u;
                        history(iter).dError = err_d;
                        history(iter).uResidual = currentUResidual;
                        history(iter).dResidual = phaseResidual;
                        history(iter).energy = totalEnergy;
                        history(iter).phaseConverged = true;
                        history(iter).dispConverged = dispConverged;

                        converged = true;
                        msg = sprintf('Converged in %d iterations (after displacement solve)', iter);
                        fprintf('Phase residual = %.6e < tolerance %.6e, and displacement converged\n', phaseResidual, d_tol);
                        fprintf('\n=== STAGGERED SOLVER CONVERGED ===\n');
                        return;
                    end

                fprintf('Step 2: Solving phase field...\n');
                [d, phaseHistory, phaseConverged, phaseMsg] = ...
                    PhaseFieldFractureSolver.solvePhaseField(mesh, mat, u, d, 3e-2, 200);

                if ~phaseConverged
                    converged = false;
                    msg = sprintf('Phase field solver failed to converge: %s', phaseMsg);
                    fprintf('\n=== STAGGERED SOLVER: Phase field solver failed ===\n');
                    return;
                end

                u_bc = PhaseFieldFractureSolver.applyBoundaryConditions(u, fixedDOFs, prescribedDOFs, prescribedValues);
                newUResidual = norm(PhaseFieldFractureSolver.assembleDisplacementResidual(mesh, mat, u_bc, d));
                finalDResidual = norm(PhaseFieldFractureSolver.assemblePhaseResidual(mesh, mat, u, d));

                err_u = norm(u - u_prev);
                err_d = norm(d - d_prev);
                totalEnergy = PhaseFieldFractureSolver.computeTotalEnergy(mesh, mat, u, d);

                history(iter).iter = iter;
                history(iter).uError = err_u;
                history(iter).dError = err_d;
                history(iter).uResidual = newUResidual;
                history(iter).dResidual = finalDResidual;
                history(iter).energy = totalEnergy;
                history(iter).phaseConverged = phaseConverged;
                history(iter).dispConverged = dispConverged;

                fprintf('||Delta u|| = %.6e, ||Delta d|| = %.6e\n', err_u, err_d);
                fprintf('Energy = %.6e\n', totalEnergy);
                fprintf('U residual = %.6e (tol: %.6e), D residual = %.6e (tol: %.6e)\n', ...
                    newUResidual, u_tol, finalDResidual, d_tol);

                if newUResidual < u_tol && phaseConverged
                    converged = true;
                    msg = sprintf('Converged in %d iterations (after phase field solve)', iter);
                    fprintf('\n=== STAGGERED SOLVER CONVERGED ===\n');
                    return;
                end

                if strcmp(dispMsg, 'NaN or Inf in increment') || ...
                   strcmp(dispMsg, 'Singular matrix') || ...
                   strcmp(phaseMsg, 'NaN or Inf in increment')
                    converged = false;
                    msg = sprintf('Numerical failure: disp=%s, phase=%s', dispMsg, phaseMsg);
                    fprintf('\n=== STAGGERED SOLVER FAILED ===\n');
                    return;
                end

                if iter == maxIterations
                    converged = false;
                    msg = sprintf('Max iterations (%d) reached', maxIterations);
                    fprintf('\n=== STAGGERED SOLVER: Max iterations reached ===\n');
                    return;
                end
            end
        end

        function u_bc = applyBoundaryConditions(u, fixedDOFs, prescribedDOFs, prescribedValues)
            % applyBoundaryConditions - Apply displacement boundary conditions
            %   u_bc = applyBoundaryConditions(u, fixedDOFs, prescribedDOFs, prescribedValues)
            %   applies boundary conditions to the displacement vector.
            %
            %   Input:
            %       u                - displacement vector
            %       fixedDOFs       - indices of DOFs with zero displacement
            %       prescribedDOFs  - indices of DOFs with prescribed displacement
            %       prescribedValues - values at prescribed DOFs
            %
            %   Output:
            %       u_bc - displacement vector with boundary conditions applied
            %
            %   Note: DOF indices are checked to ensure they are within bounds
            u_bc = u;
            validFixed = fixedDOFs(fixedDOFs <= length(u));
            u_bc(validFixed) = 0;
            validPrescribed = prescribedDOFs(prescribedDOFs <= length(u));
            u_bc(validPrescribed) = prescribedValues(1:length(validPrescribed));
        end
    end

    %===========================================================================
    %                        POST PROCESSING
    %===========================================================================
    methods (Static)
        function plotDisplacement(mesh, u, style)
            % plotDisplacement - Plot displacement field
            %   plotDisplacement(mesh, u) creates a contour plot of the displacement field.
            %
            %   Input:
            %       mesh - mesh structure
            %       u    - displacement vector
            %       style - plotting style: 'contourf' (default) or 'scatter'
            %
            %   See also: plotPhaseField, plotResults
            if nargin < 3
                style = 'contourf';
            end

            figure('Name', 'Displacement Field', 'NumberTitle', 'off');
            if strcmp(style, 'scatter')
                scatter3(mesh.x(:), mesh.y(:), u(:), 10, u(:));
                xlabel('X'); ylabel('Y'); zlabel('Displacement');
                colorbar;
            else
                contourf(mesh.x, mesh.y, reshape(u, mesh.nx, mesh.ny));
                xlabel('X'); ylabel('Y');
                colorbar;
            end
            title('Displacement Field');
        end

        function plotPhaseField(mesh, d, style)
            % plotPhaseField - Plot phase field (damage) distribution
            %   plotPhaseField(mesh, d) creates a contour plot of the phase field.
            %
            %   Input:
            %       mesh - mesh structure
            %       d    - phase field vector (0 = intact, 1 = fully broken)
            %       style - plotting style: 'contourf' (default) or 'scatter'
            %
            %   See also: plotDisplacement, plotResults
            if nargin < 3
                style = 'contourf';
            end

            figure('Name', 'Phase Field', 'NumberTitle', 'off');
            if strcmp(style, 'scatter')
                scatter3(mesh.x(:), mesh.y(:), d(:), 10, d(:));
                xlabel('X'); ylabel('Y'); zlabel('Damage');
                colorbar;
            else
                contourf(mesh.x, mesh.y, reshape(d, mesh.nx, mesh.ny));
                xlabel('X'); ylabel('Y');
                colorbar;
            end
            title('Phase Field (Damage)');
        end

        function plotStress(mesh, mat, u, d)
            % plotStress - Plot stress field
            %   plotStress(mesh, mat, u, d) computes and plots the stress field
            %   based on the displacement and phase field solutions.
            %
            %   Input:
            %       mesh - mesh structure
            %       mat  - material properties
            %       u    - displacement vector
            %       d    - phase field vector
            %
            %   See also: plotDisplacement, plotPhaseField
            stress = zeros(mesh.nx-1, mesh.ny-1);

            for e = 1:mesh.nElements
                nodes = mesh.elements(e, :);
                elem_i = mod(double(e)-1, double(mesh.nx)-1) + 1;
                elem_j = floor((double(e)-1) / (double(mesh.nx)-1)) + 1;

                xi = 0; eta = 0;
                N = PhaseFieldFractureSolver.shapeFunctions(xi, eta);
                [dN_dxi, dN_deta] = PhaseFieldFractureSolver.shapeFunctionDerivatives(xi, eta);
                [~, ~, ~, dN_dx, dN_dy] = PhaseFieldFractureSolver.computeJacobian(mesh, nodes, dN_dxi, dN_deta);
                B = PhaseFieldFractureSolver.computeBMatrix(dN_dx, dN_dy);

                u_elem = u(nodes);
                strain = B * u_elem;
                d_val = dot(N, d(nodes));
                g = PhaseFieldFractureSolver.degradationFunction(d_val);

                sigma = g * mat.G * strain;
                stress(elem_i, elem_j) = sqrt(sigma(1)^2 + sigma(2)^2);
            end

            figure('Name', 'Shear Stress', 'NumberTitle', 'off');
            contourf(mesh.x(1:end-1, 1:end-1)', mesh.y(1:end-1, 1:end-1)', stress');
            xlabel('X'); ylabel('Y');
            colorbar;
            title('Shear Stress Magnitude');
        end

        function plotHistory(loadSteps, displacements, damages, energies)
            figure('Name', 'Convergence History', 'NumberTitle', 'off');

            subplot(2, 2, 1);
            plot(loadSteps, displacements, '-o');
            xlabel('Load Step'); ylabel('Max Displacement');
            title('Displacement History');

            subplot(2, 2, 2);
            plot(loadSteps, damages, '-o');
            xlabel('Load Step'); ylabel('Max Damage');
            title('Damage History');

            subplot(2, 2, 3);
            plot(loadSteps, energies, '-o');
            xlabel('Load Step'); ylabel('Total Energy');
            title('Energy History');

            subplot(2, 2, 4);
            plot(displacements, damages, '-o');
            xlabel('Displacement'); ylabel('Damage');
            title('Force-Damage Curve');
        end

        function saveToVTK(mesh, mat, u, d, filename)
            fid = fopen([filename, '.vtk'], 'w');
            fprintf(fid, '# vtk DataFile Version 3.0\n');
            fprintf(fid, 'Phase Field Fracture Results\n');
            fprintf(fid, 'ASCII\n');
            fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
            fprintf(fid, 'DIMENSIONS %d %d 1\n', mesh.nx, mesh.ny);
            fprintf(fid, 'ORIGIN 0 0 0\n');
            fprintf(fid, 'SPACING %.6f %.6f 1\n', mesh.dx, mesh.dy);
            fprintf(fid, 'POINT_DATA %d\n', mesh.nNodes);

            fprintf(fid, 'SCALARS displacement double\n');
            fprintf(fid, 'LOOKUP_TABLE default\n');
            for j = 1:mesh.ny
                for i = 1:mesh.nx
                    fprintf(fid, '%.6e\n', u((j-1)*mesh.nx + i));
                end
            end

            fprintf(fid, 'SCALARS phase_field double\n');
            fprintf(fid, 'LOOKUP_TABLE default\n');
            for j = 1:mesh.ny
                for i = 1:mesh.nx
                    fprintf(fid, '%.6e\n', d((j-1)*mesh.nx + i));
                end
            end

            fprintf(fid, 'SCALARS history_variable double\n');
            fprintf(fid, 'LOOKUP_TABLE default\n');
            for j = 1:mesh.ny
                for i = 1:mesh.nx
                    h_val = 0;
                    elem_i_min = max(1, i-1);
                    elem_i_max = min(mesh.nx-1, i);
                    elem_j_min = max(1, j-1);
                    elem_j_max = min(mesh.ny-1, j);
                    for ei = elem_i_min:elem_i_max
                        for ej = elem_j_min:elem_j_max
                            h_val = max(h_val, max(mesh.H(ei, ej, :)));
                        end
                    end
                    fprintf(fid, '%.6e\n', h_val);
                end
            end

            fprintf(fid, 'SCALARS history_variable_prev double\n');
            fprintf(fid, 'LOOKUP_TABLE default\n');
            for j = 1:mesh.ny
                for i = 1:mesh.nx
                    h_val = 0;
                    elem_i_min = max(1, i-1);
                    elem_i_max = min(mesh.nx-1, i);
                    elem_j_min = max(1, j-1);
                    elem_j_max = min(mesh.ny-1, j);
                    for ei = elem_i_min:elem_i_max
                        for ej = elem_j_min:elem_j_max
                            h_val = max(h_val, max(mesh.Hprev(ei, ej, :)));
                        end
                    end
                    fprintf(fid, '%.6e\n', h_val);
                end
            end

            % 计算并输出应力分量（反平面剪切问题：τ_xz 和 τ_yz）
            tau_xz = zeros(mesh.nx, mesh.ny);
            tau_yz = zeros(mesh.nx, mesh.ny);
            
            for e = 1:mesh.nElements
                nodes = mesh.elements(e, :);
                elem_i = mod(double(e)-1, double(mesh.nx)-1) + 1;
                elem_j = floor((double(e)-1) / (double(mesh.nx)-1)) + 1;

                xi = 0; eta = 0;
                N = PhaseFieldFractureSolver.shapeFunctions(xi, eta);
                [dN_dxi, dN_deta] = PhaseFieldFractureSolver.shapeFunctionDerivatives(xi, eta);
                [~, ~, ~, dN_dx, dN_dy] = PhaseFieldFractureSolver.computeJacobian(mesh, nodes, dN_dxi, dN_deta);
                B = PhaseFieldFractureSolver.computeBMatrix(dN_dx, dN_dy);

                u_elem = u(nodes);
                strain = B * u_elem;
                d_val = dot(N, d(nodes));
                g = PhaseFieldFractureSolver.degradationFunction(d_val);

                sigma = g * mat.G * strain;
                s_xz = sigma(1);
                s_yz = sigma(2);
                
                tau_xz(elem_i:elem_i+1, elem_j:elem_j+1) = max(tau_xz(elem_i:elem_i+1, elem_j:elem_j+1), abs(s_xz));
                tau_yz(elem_i:elem_i+1, elem_j:elem_j+1) = max(tau_yz(elem_i:elem_i+1, elem_j:elem_j+1), abs(s_yz));
            end

            fprintf(fid, 'SCALARS shear_stress_xz double\n');
            fprintf(fid, 'LOOKUP_TABLE default\n');
            for j = 1:mesh.ny
                for i = 1:mesh.nx
                    fprintf(fid, '%.6e\n', tau_xz(i, j));
                end
            end

            fprintf(fid, 'SCALARS shear_stress_yz double\n');
            fprintf(fid, 'LOOKUP_TABLE default\n');
            for j = 1:mesh.ny
                for i = 1:mesh.nx
                    fprintf(fid, '%.6e\n', tau_yz(i, j));
                end
            end

            % 计算并输出等效应力（von Mises）
            stress_magnitude = sqrt(tau_xz.^2 + tau_yz.^2);
            fprintf(fid, 'SCALARS von_mises_stress double\n');
            fprintf(fid, 'LOOKUP_TABLE default\n');
            for j = 1:mesh.ny
                for i = 1:mesh.nx
                    fprintf(fid, '%.6e\n', stress_magnitude(i, j));
                end
            end

            fclose(fid);
            fprintf('VTK file saved: %s.vtk\n', filename);
        end

        function saveResultsToText(mesh, u, d, filename)
            fid = fopen([filename, '_displacement.txt'], 'w');
            for j = 1:mesh.ny
                for i = 1:mesh.nx
                    fprintf(fid, '%.6e ', u((j-1)*mesh.nx + i));
                end
                fprintf(fid, '\n');
            end
            fclose(fid);
            fprintf('Displacement saved: %s_displacement.txt\n', filename);

            fid = fopen([filename, '_phase.txt'], 'w');
            for j = 1:mesh.ny
                for i = 1:mesh.nx
                    fprintf(fid, '%.6e ', d((j-1)*mesh.nx + i));
                end
                fprintf(fid, '\n');
            end
            fclose(fid);
            fprintf('Phase field saved: %s_phase.txt\n', filename);
        end
    end
end