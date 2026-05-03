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

function runRectangleFractureExample()
    % RECTANGLEFRACTUREEXAMPLE 长方形边值问题算例
    % 用于模拟PMMA材料在反平面剪切载荷下的断裂行为
    
    fprintf('========================================\n');
    fprintf('Rectangle Boundary Value Problem Example\n');
    fprintf('Phase Field Fracture Simulation\n');
    fprintf('========================================\n\n');

    %===========================================================================
    % 1. 几何模型定义
    %===========================================================================
    fprintf('1. 几何模型定义\n');
    fprintf('-----------------\n');
    
    Lx = 2.0;    % 长方形区域长度 (m)
    Ly = 1.0;    % 长方形区域高度 (m)
    nx = 61;     % x方向节点数 (dx = 0.033 m, 约1800单元)
    ny = 31;     % y方向节点数 (dy = 0.033 m)
    
    fprintf('  区域尺寸: Lx = %.2f m, Ly = %.2f m\n', Lx, Ly);
    fprintf('  网格划分: nx = %d, ny = %d\n', nx, ny);
    fprintf('  单元尺寸: dx = %.6f m, dy = %.6f m\n', Lx/(nx-1), Ly/(ny-1));
    fprintf('  总节点数: %d\n', nx * ny);
    fprintf('  总单元数: %d\n', (nx-1) * (ny-1));

    % 创建网格
    mesh = PhaseFieldFractureSolver.createMesh(nx, ny, Lx, Ly);

    %===========================================================================
    % 2. 材料属性 (PMMA)
    %===========================================================================
    fprintf('\n2. 材料属性 (PMMA)\n');
    fprintf('-------------------\n');
    
    E = 3.0e9;   % 弹性模量 (Pa) - PMMA典型值
    nu = 0.35;   % 泊松比 - PMMA典型值
    Gc = 450;    % 临界能量释放率 (J/m^2) - PMMA典型值
    l0 = 0.08;   % 长度尺度参数 (m) - 满足 l0 >= 2*dx (0.08 >= 0.067)
    
    fprintf('  弹性模量 E = %.2e Pa\n', E);
    fprintf('  泊松比 nu = %.2f\n', nu);
    fprintf('  剪切模量 G = %.2e Pa\n', E/(2*(1+nu)));
    fprintf('  临界能量释放率 Gc = %.2e J/m^2\n', Gc);
    fprintf('  长度尺度参数 l0 = %.4f m\n', l0);
    
    % 创建材料
    mat = PhaseFieldFractureSolver.createMaterial(E, nu, Gc, l0);

    %===========================================================================
    % 3. 初始条件 - 预制裂纹
    %===========================================================================
    fprintf('\n3. 初始条件 - 预制裂纹\n');
    fprintf('-----------------------\n');
    
    d_prev = zeros(double(mesh.nPhase), 1);
    
    % 定义裂纹区域：中间高度位置，左侧非零值，右侧零值
    crack_center_y = Ly / 2;       % 裂纹中心高度
    crack_width = 0.02;            % 裂纹宽度 (m)
    crack_length = Lx / 2;         % 裂纹长度 - 从左侧到中间
    
    fprintf('  裂纹中心高度: y = %.2f m\n', crack_center_y);
    fprintf('  裂纹宽度: %.4f m\n', crack_width);
    fprintf('  裂纹长度: %.2f m (从左侧边界到中间位置)\n', crack_length);
    fprintf('  初始损伤值: d = 0.95 (裂纹区域)\n');
    
    % 设置初始相场分布
    for i = 1:mesh.nx
        for j = 1:mesh.ny
            x = mesh.x(i, j);
            y = mesh.y(i, j);
            
            % 判断是否在裂纹区域内
            in_crack_region = (x <= crack_length) && ...  % 左侧区域
                             (abs(y - crack_center_y) <= crack_width/2);
            
            if in_crack_region
                n = (j-1)*mesh.nx + i;
                d_prev(n) = 0.95;  % 高损伤值表示预制裂纹
            end
        end
    end
    
    % 显示初始裂纹信息
    fprintf('  初始相场统计: min(d) = %.4f, max(d) = %.4f\n', min(d_prev), max(d_prev));
    fprintf('  裂纹区域节点数: %d\n', sum(d_prev > 0.5));
    
    % 设置初始历史变量H（在裂纹区域设置非零值）
    initial_H_value = Gc / l0;  % 初始历史变量值 = Gc/l0（临界能量释放率除以长度尺度参数）
    fprintf('  初始H值设置为 Gc/l0 = %.2e / %.4f = %.2e\n', Gc, l0, initial_H_value);
    
    crack_element_count = 0;
    for elem_i = 1:mesh.nx-1
        for elem_j = 1:mesh.ny-1
            % 单元中心坐标 = 单元左下角节点坐标 + 半网格间距
            x_center = mesh.x(elem_i, elem_j) + mesh.dx/2;
            y_center = mesh.y(elem_i, elem_j) + mesh.dy/2;
            
            % 判断是否在裂纹区域内（从左到中间，中间高度附近）
            % 使用网格间距作为裂纹宽度的下限，确保至少有一个单元被选中
            effective_crack_width = max(crack_width, mesh.dy);
            in_crack_region = (x_center <= crack_length) && ...
                             (abs(y_center - crack_center_y) <= effective_crack_width/2);
            
            if in_crack_region
                mesh.H(elem_i, elem_j, :) = initial_H_value;
                mesh.Hprev(elem_i, elem_j, :) = initial_H_value;
                crack_element_count = crack_element_count + 1;
            end
        end
    end
    fprintf('  裂纹区域单元数: %d\n', crack_element_count);
    fprintf('  初始历史变量H: min=%.6e, max=%.6e\n', min(mesh.H(:)), max(mesh.H(:)));
    fprintf('  初始历史变量Hprev: min=%.6e, max=%.6e\n', min(mesh.Hprev(:)), max(mesh.Hprev(:)));

    %===========================================================================
    % 4. 边界条件
    %===========================================================================
    fprintf('\n4. 边界条件\n');
    fprintf('-----------\n');
    
    u_prev = zeros(double(mesh.nDOF), 1);
    % d_prev 已经在前面设置了预制裂纹，这里不再重新初始化
    
    % 保存初始状态到VTK文件
    vtkFilenameInitial = 'rectangle_fracture_initial';
    PhaseFieldFractureSolver.saveToVTK(mesh, mat, u_prev, d_prev, vtkFilenameInitial);
    fprintf('  初始状态VTK文件已保存: %s.vtk\n', vtkFilenameInitial);
    
    fixedDOFs = [];
    prescribedDOFs = [];
    prescribedValues = [];
    
    % 反平面剪切边界条件:
    % 左侧边界: 自由
    % 右侧边界: 自由
    % 上部边界: 施加正向位移约束 w = +U
    % 下部边界: 施加负向位移约束 w = -U
    % 单点约束: 消除刚体位移模式
    
    fprintf('  左侧边界 (x = %.2f m): 自由\n', 0);
    fprintf('  右侧边界 (x = %.2f m): 自由\n', Lx);
    fprintf('  上部边界 (y = %.2f m): 正向位移约束 w = +U\n', Ly);
    fprintf('  下部边界 (y = %.2f m): 负向位移约束 w = -U\n', 0);
    fprintf('  单点约束: 中心点固定 w = 0 (消除刚体位移)\n');
    
    % 单点约束 - 中心点固定，消除刚体位移模式
    center_node = round(mesh.ny/2)*mesh.nx + round(mesh.nx/2);
    fixedDOFs = [fixedDOFs, center_node];
    
    % 上部边界节点 (y = Ly)
    for i = 1:mesh.nx
        n = (mesh.ny-1)*mesh.nx + i;
        prescribedDOFs = [prescribedDOFs, n];
        prescribedValues = [prescribedValues, 0];  % 初始值
    end
    
    % 下部边界节点 (y = 0)
    for i = 1:mesh.nx
        n = i;
        prescribedDOFs = [prescribedDOFs, n];
        prescribedValues = [prescribedValues, 0];  % 初始值
    end
    
    fprintf('  固定节点数: %d\n', length(fixedDOFs));
    fprintf('  位移约束节点数: %d\n', length(prescribedDOFs));

    %===========================================================================
    % 5. 加载参数
    %===========================================================================
    fprintf('\n5. 加载参数\n');
    fprintf('-----------\n');
    
    % 调整前参数: U_max = 1e-3 m, deltaU = 1e-5 m
    % 调整后参数: 基于收敛分析，将最大位移减小到能够稳定收敛的范围
    U_max = 5e-5;        % 最大位移载荷 (m) - 调整为能收敛的最大值
    deltaU = 1e-6;       % 位移增量 (m) - 减小步长保持加载稳定性
    nSteps = ceil(U_max / deltaU);
    
    fprintf('  目标位移: U_max = %.2e m (调整前: 1.00e-03 m)\n', U_max);
    fprintf('  位移增量: deltaU = %.2e m (调整前: 1.00e-05 m)\n', deltaU);
    fprintf('  加载步数: nSteps = %d\n', nSteps);

    %===========================================================================
    % 6. 求解器设置
    %===========================================================================
    fprintf('\n6. 求解器设置\n');
    fprintf('-------------\n');
    
    tolerance = 1e-6;        % 收敛容差（调整为1e-6）
    maxIterations = 100;     % 最大交错迭代次数
    
    fprintf('  收敛容差: %.1e\n', tolerance);
    fprintf('  最大交错迭代次数: %d\n', maxIterations);
    fprintf('  求解器类型: 交错Newton-Raphson（无线搜索）\n');
    fprintf('  相场模型: AT2 (Ambati et al. 2015)\n');

    %===========================================================================
    % 7. 结果存储
    %===========================================================================
    loadSteps = zeros(nSteps, 1);
    displacements = zeros(nSteps, 1);
    damages = zeros(nSteps, 1);
    energies = zeros(nSteps, 1);
    H_max_history = zeros(nSteps, 1);

    %===========================================================================
    % 8. 主求解循环
    %===========================================================================
    fprintf('\n7. 开始求解\n');
    fprintf('-----------\n');
    
    tic;
    
    for step = 1:nSteps
        U = step * deltaU;
        
        % 更新加载步
        mesh = PhaseFieldFractureSolver.startNewLoadStep(mesh);
        
        fprintf('\n========== Load Step %d/%d: U = %.2e ==========\n', step, nSteps, U);
        
        u0 = u_prev;
        d0 = d_prev;
        
        % 更新位移边界条件
        % 上部边界: 正向位移 +U
        % 下部边界: 负向位移 -U
        stepPrescribedValues = zeros(size(prescribedDOFs));
        % 上部边界节点 (前mesh.nx个prescribedDOFs)
        stepPrescribedValues(1:mesh.nx) = U;      % 上部: +U
        % 下部边界节点 (后mesh.nx个prescribedDOFs)
        stepPrescribedValues(mesh.nx+1:end) = -U;  % 下部: -U
        
        % 调用交错求解器
        [u, d, mesh, history, converged, msg] = ...
            PhaseFieldFractureSolver.solveStaggered(...
                mesh, mat, u0, d0, ...
                fixedDOFs, prescribedDOFs, stepPrescribedValues, ...
                tolerance, maxIterations);
        
        if ~converged
            fprintf('\n========== SIMULATION TERMINATED ==========\n');
            fprintf('Convergence failure at step %d\n', step);
            fprintf('Message: %s\n', msg);
            break;
        end
        
        % 保存当前加载步的结果到VTK文件
        vtkFilenameStep = sprintf('rectangle_fracture_step_%03d', step);
        PhaseFieldFractureSolver.saveToVTK(mesh, mat, u, d, vtkFilenameStep);
        fprintf('  VTK file saved for step %d: %s.vtk\n', step, vtkFilenameStep);
        
        % 更新状态
        u_prev = u;
        d_prev = d;
        
        % 存储结果
        loadSteps(step) = U;
        displacements(step) = max(abs(u));
        damages(step) = max(d);
        energies(step) = history(end).energy;
        H_max_history(step) = max(mesh.H(:));
        
        % 输出当前步结果
        fprintf('  最大位移: %.6e m\n', displacements(step));
        fprintf('  最大损伤: %.6f\n', damages(step));
        fprintf('  总能量: %.6e J\n', energies(step));
        fprintf('  历史变量H: min=%.6e, max=%.6e\n', min(mesh.H(:)), max(mesh.H(:)));
        fprintf('  收敛迭代次数: %d\n', length([history.iter]));
    end
    
    elapsedTime = toc;
    
    %===========================================================================
    % 9. 结果输出
    %===========================================================================
    fprintf('\n8. 结果输出\n');
    fprintf('-----------\n');
    
    fprintf('\n========================================\n');
    fprintf('Simulation Complete\n');
    fprintf('========================================\n');
    fprintf('CPU Time: %.2f seconds\n', elapsedTime);
    
    if ~isempty(loadSteps)
        fprintf('\nFinal Results:\n');
        fprintf('  Final displacement: %.6e m\n', displacements(end));
        fprintf('  Final max damage: %.6f\n', damages(end));
        fprintf('  Final energy: %.6e J\n', energies(end));
    end

    %===========================================================================
    % 10. 可视化输出
    %===========================================================================
    fprintf('\n9. 可视化输出\n');
    fprintf('-------------\n');
    
    % 保存VTK文件到MATLAB目录
    vtkFilename = 'rectangle_fracture_simulation';
    PhaseFieldFractureSolver.saveToVTK(mesh, mat, u, d, vtkFilename);
    fprintf('  VTK文件已保存: %s.vtk\n', vtkFilename);
    
    % 保存文本结果到MATLAB目录
    txtFilename = 'rectangle_fracture_results';
    PhaseFieldFractureSolver.saveResultsToText(mesh, u, d, txtFilename);
    fprintf('  文本结果已保存: %s_displacement.txt, %s_phase.txt\n', txtFilename, txtFilename);
    
    % 绘制收敛曲线
    figure('Name', 'Convergence History', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 400]);
    subplot(1, 4, 1);
    plot(loadSteps(1:step), displacements(1:step), 'b-o', 'LineWidth', 2);
    xlabel('Displacement (m)'); ylabel('Max Displacement (m)');
    title('Displacement vs Load');
    grid on;
    
    subplot(1, 4, 2);
    plot(loadSteps(1:step), damages(1:step), 'r-o', 'LineWidth', 2);
    xlabel('Displacement (m)'); ylabel('Max Damage');
    title('Damage vs Load');
    grid on;
    
    subplot(1, 4, 3);
    plot(loadSteps(1:step), energies(1:step), 'g-o', 'LineWidth', 2);
    xlabel('Displacement (m)'); ylabel('Energy (J)');
    title('Energy vs Load');
    grid on;
    
    subplot(1, 4, 4);
    plot(loadSteps(1:step), H_max_history(1:step), 'm-o', 'LineWidth', 2);
    xlabel('Displacement (m)'); ylabel('Max History Variable');
    title('History Variable vs Load');
    grid on;
    
    saveas(gcf, 'convergence_history.png');
    fprintf('  收敛曲线已保存: convergence_history.png\n');
    
    % 绘制最终结果云图
    figure('Name', 'Final Results', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 500]);
    
    subplot(1, 3, 1);
    PhaseFieldFractureSolver.plotDisplacement(mesh, u, 'scatter');
    title('Displacement Field (w)');
    
    subplot(1, 3, 2);
    PhaseFieldFractureSolver.plotPhaseField(mesh, d, 'scatter');
    title('Phase Field (Damage)');
    
    subplot(1, 3, 3);
    PhaseFieldFractureSolver.plotStress(mesh, mat, u, d);
    title('Shear Stress (τ_xy)');
    
    saveas(gcf, 'final_results.png');
    fprintf('  结果云图已保存: final_results.png\n');
    
    fprintf('\n所有结果已保存至 MATLAB 目录\n');
end
