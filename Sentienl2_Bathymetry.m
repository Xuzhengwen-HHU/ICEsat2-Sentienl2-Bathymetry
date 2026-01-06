clear
clc
close all

%% 读取Icesat2水深数据集
cd L:\Icesat2_xzw\Cor_data
filename1='Icesat2_corrected_merged5.csv';
% filename1='Icesat2_corrected_20230506.csv';
% filename1='Icesat2_corrected_5789.csv';
% filename1='Icesat2_corrected_1456789.csv';

Data=importdata(filename1).data;
Depth=[Data(:,2),Data(:,1),Data(:,5)];%% 经度，纬度，水深
%%%做反演时刻的潮汐校正
% %%Sentienl2日期潮汐（世界时To北京时+8h,超过24h ,天数+1）
tide_20160217_11=(0.78+0.85)/2-0.95;
tide_20190224_11=(0.69+0.81)/2-0.95;
tide_20220315_11=(0.86+0.92)/2-0.95;
tide_20240208_11=(0.73+0.83)/2-0.95;

Depth=[Data(:,2),Data(:,1),Data(:,5)-tide_20190224_11];
data=Depth;%生成水深点位置及水深信息矩阵, 经度，纬度，水深
%% 读取反演时刻的反射率
% filename2='L:\Sentienl2_xzw\fanshelv\Sen2_20160217_fanshelv.txt';
filename2='L:\Sentienl2_xzw\fanshelv\Sen2_20190224_fanshelv.txt';
% filename2='L:\Sentienl2_xzw\fanshelv\Sen2_20220315_fanshelv.txt';
% filename2='L:\Sentienl2_xzw\fanshelv\Sen2_20240208_fanshelv.txt';

fanshelv=importdata(filename2);

fanshelv=fanshelv.data;%% 纬度，经度，反射率
dataname=filename2(32:39);
dataname=(['Data : ', dataname(1:4),'-',dataname(5:6),'-',dataname(7:8) ] );
%% 读取验证水深数据
filename3 = 'L:\code\Si_tu\si_tu_depth.txt';
Data_situ=importdata(filename3).data;%%%%整个永乐环礁的反演数据
Data_situ=[Data_situ(17061:20552,:);Data_situ(26387:30282,:)];

lon = Data_situ(:,1);
lat = Data_situ(:,2);
hei = Data_situ(:,3);
%%%做反演时刻的潮汐校正
hei = hei-tide_20190224_11;

situ_data=[lon,lat,hei];

%% 水深点重采样（匹配Sentienl2和Icesat2数据）
%%%%由于Sentienl2为10米分辨率的栅格图像，一个栅格中包含多个Icesat2光子点
rows=2243;cols=3654;
lat1=min(fanshelv(:,1));lat2=max(fanshelv(:,1));lat_dt=(lat2-lat1)/rows;
lon1=min(fanshelv(:,2));lon2=max(fanshelv(:,2));lon_dt=(lon2-lon1)/cols;
latitude = fanshelv(:, 1); % 第一列为纬度
longitude = fanshelv(:, 2); % 第二列为经度

[XI,YI] = meshgrid(lon1:lon_dt:lon2,lat1:lat_dt:lat2);%构建网格
x_num=size(XI,1);%读取划分网格的横坐标数2244
y_num=size(XI,2);%读取划分网格的纵坐标数3655
for x=1:y_num-1%length(XI(1,:))%重采样水深横坐标，依据22个标记点，划分为21间隔循环切片提取
    ll=find (data(:,1)>XI(1,x)&data(:,1)<XI(1,x+1));
   for y=1:x_num-1%length(YI(:,1)) %重采样水深纵坐标，依据243个标记点，划分为242间隔循环切片提取
         hh=find(data(ll,2)>YI(y,1) &data(ll,2)<YI(y+1,1));
         T(y,x)=mean(data(ll(hh),3));%以一列列数据存储重采样水深值
         Lat_resample(y,x)=mean(data(ll(hh),2));%以一列列数据存储重采样纬度值
   end
end
T=flipud(T);
imshow(T,[])
%% 验证水深点重采样
for x=1:y_num-1%length(XI(1,:))%重采样水深横坐标，依据22个标记点，划分为21间隔循环切片提取
    ll=find (situ_data(:,1)>XI(1,x)&situ_data(:,1)<XI(1,x+1));
   for y=1:x_num-1%length(YI(:,1)) %重采样水深纵坐标，依据243个标记点，划分为242间隔循环切片提取
         hh=find(situ_data(ll,2)>YI(y,1) &situ_data(ll,2)<YI(y+1,1));
         situ_Depth(y,x)=mean(situ_data(ll(hh),3));%以一列列数据存储重采样水深值
         Lat_resample(y,x)=mean(situ_data(ll(hh),2));%以一列列数据存储重采样纬度值
   end
end
situ_Depth=flipud(situ_Depth);
% situ_Depth=situ_Depth(situ_Depth<20);
% imshow(situ_Depth,[])

%% stumpf比值波段模型建立(用Icesat2所测水深与Sentienl2反射率构建模型，得到反演参数)
if length(fanshelv)~=3654*2243;
    cols=cols+1;
end

Band_blue=(reshape(fanshelv(:,3),cols,rows))';
Band_green=(reshape(fanshelv(:,4),cols,rows))';
Band_red=(reshape(fanshelv(:,5),cols,rows))';
Band_nir=(reshape(fanshelv(:,6),cols,rows))';
if length(fanshelv)~=3654*2243
    Band_blue=Band_blue(:,2:end);
    Band_green=Band_green(:,2:end);
    Band_red=Band_red(:,2:end);
    Band_nir=Band_nir(:,2:end);
end

ln_b=log(Band_blue);
% imshow(ln_b,[])
ln_g=log(Band_green);
ln_r=log(Band_red);
ln_nir=log(Band_nir);
% imshow(Band_blue,[])
% 合并成RGB图像
rgbImage = cat(3, Band_red, Band_green, Band_blue);
% 将数据归一化到0-255范围（如果数据值不在0-255范围内）
rgbImage1 = uint8(mat2gray(rgbImage) * 255);
rgbImage=imadjust(rgbImage1,[0 0 0;0.2 0.35 0.35],[]);  %imadjust()对RGB图像进行处理
% %%%%%%%%%%永乐环礁卫星图
% figure(3)
% % set(gcf,'position',[400 250 600 380])
% % set(gca,'position',[0.03 0.05 0.9 0.95]);
% % 设置地图投影和显示范围
% m_proj('Mercator', 'longitudes', [min(longitude) max(longitude)], 'latitudes', [min(latitude) max(latitude)]);
% % 插入遥感图（注意纬度顺序要颠倒）
% m_image([min(longitude) max(longitude)], [max(latitude) min(latitude)], rgbImage);
% hold on
% % 自定义经纬度刻度
% xticks = 111.5:0.05:111.8;
% yticks = 16.45:0.05:16.60;
% 
% xtick_labels = arrayfun(@(x) sprintf('%.2f°E', x), xticks, 'UniformOutput', false);
% ytick_labels = arrayfun(@(y) sprintf('%.2f°N', y), yticks, 'UniformOutput', false);
% % 添加网格线和刻度标签
% m_grid('box', 'on', 'tickdir', 'out', ...
%        'XaxisLocation', 'bottom', 'YaxisLocation', 'right', ...
%        'fontsize', 8, 'fontname', 'Times New Roman', ...
%        'tickstyle', 'dd', ...
%        'xtick', xticks, 'ytick', yticks, ...
%        'xticklabels', xtick_labels, 'yticklabels', ytick_labels);
% % 在第 1750 行对应的纬度处画红色水平线
% latitude_matrix=(reshape(latitude,cols,rows))';
% lat_cut = latitude_matrix(1750);
% m_line([min(longitude), max(longitude)], ...   % 起点和终点经度
%        [lat_cut, lat_cut], ...                   % 相同纬度
%        'Color','r','LineWidth',1.5);
% 
% 
% 
% 
% 
% x_blue=1:3654;
% % Prepare longitude vector
% n = numel(x_blue);
% lon = linspace(min(longitude) , max(longitude), n);
% 
% figure;
% set(gcf, 'Color', 'w');
% hold on;
% % --- Left Y-axis: setup and shading ---
% yyaxis left
% ylim([0 0.55]);       % reflectance limits
% yl = ylim;
% % shading between 111.596–111.628
% patch([111.57, 111.607, 111.607, 111.57], ...
%       [yl(1),    yl(1),   0.45,    0.45], ...
%       [0.6 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% % shading between 111.743–111.766
% patch([111.736, 111.763, 111.763, 111.736], ...
%       [yl(1),    yl(1), 0.45,0.45], ...
%       [0.6 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% % plot reflectance curves
% h1 = plot(lon, Band_blue(1750,:)./10000,  '-b', 'LineWidth', 1.5);
% h2 = plot(lon, Band_green(1750,:)./10000, '-g', 'LineWidth', 1.5);
% ylabel('Reflectance', 'FontName', 'Times New Roman', 'FontSize', 14);
% ax = gca;
% ax.YAxis(1).Color = [0 0 1];
% 
% % --- Right Y-axis: log ratio ---
% yyaxis right
% xlim([min(longitude) , max(longitude)])
% ylim([0.5 1.3]);
% h3 = plot(lon, ln_b(1750,:)./ln_g(1750,:), '-k', 'LineWidth', 1.5);
% ylabel('Log Ratio (ln(B)/ln(G))', 'FontName', 'Times New Roman', 'FontSize', 14);
% ax.YAxis(2).Color = [0 0 0];
% 
% % --- X-axis formatting ---
% % xlabel('Longitude (°E)', 'FontName', 'Times New Roman', 'FontSize', 14);
% xticks = 111.5:0.1:111.8;
% set(ax, 'XTick', xticks, ...
%        'XTickLabel', arrayfun(@(x) sprintf('%.2f°E', x), xticks, 'UniformOutput', false));
% 
% % --- Legend, grid, styling ---
% legend([h1, h2, h3], {'Blue Reflectance','Green Reflectance','Log Ratio'}, ...
%        'Location', 'northeast', 'FontName', 'Times New Roman', 'FontSize', 14);
% grid on;
% ax.Box = 'on';
% ax.FontName = 'Times New Roman';
% ax.FontSize = 14;
% 
% text(111.48, 0.75, 'Deep Water', 'fontsize', 14, 'color', 'k', ...
%     'fontname', 'Times New Roman', 'FontWeight', 'bold');hold on
% text(111.66, 0.75, 'Lagon', 'fontsize', 14, 'color', 'k', ...
%     'fontname', 'Times New Roman', 'FontWeight', 'bold');
% 

% 合并遥感图和光谱曲线为 1×2 子图，并添加 (A)/(B) 标注
figure('Color','w', 'Position', [100, 100, 1200, 450]);

%% 子图 A：永乐环礁卫星图
set( subplot(1,2,1), 'Position',[0.05, 0.03, 0.43, 0.95]);

% 设置地图投影和显示范围
m_proj('Mercator', ...
       'longitudes', [min(longitude) max(longitude)], ...
       'latitudes',  [min(latitude)  max(latitude)]);
% 插入遥感图（注意纬度顺序要颠倒）
m_image([min(longitude) max(longitude)], [max(latitude) min(latitude)], rgbImage);
hold on;
% 自定义经纬度刻度
xticks = 111.5:0.1:111.8;
yticks = 16.5:0.1:16.60;
xtick_labels = arrayfun(@(x) sprintf('%.1f°E', x), xticks, 'UniformOutput', false);
ytick_labels = arrayfun(@(y) sprintf('%.1f°N', y), yticks, 'UniformOutput', false);
m_grid('box','on','tickdir','out', ...
       'XaxisLocation','bottom','YaxisLocation','left', ...
       'fontsize',12,'fontname','Times New Roman', ...
       'tickstyle','dd', ...
       'xtick',xticks,'ytick',yticks, ...
       'xticklabels',xtick_labels,'yticklabels',ytick_labels);

% 画红色水平参考线
cols = size(rgbImage,2); rows = size(rgbImage,1);
latitude_matrix = reshape(latitude, cols, rows)';  
lat_cut = latitude_matrix(1750);
m_line([min(longitude), max(longitude)], [lat_cut, lat_cut], ...
       'Color','r','LineWidth',1.5);

% 添加子图标签 (A)
m_text(111.47, 16.61,'(A)', 'FontSize', 20, 'Color', 'k', ...
     'FontWeight', 'bold', 'BackgroundColor', 'w', ...
     'Margin', 1, 'EdgeColor', 'none');


% 在各国家大致位置标注岛屿名称
m_text(111.5, 16.467, 'Jinyin Island', 'fontsize', 8, 'color', 'w', ...
    'fontname', 'Times New Roman', 'FontWeight', 'bold');
m_text(111.565, 16.435, 'Lingyang Jiao', 'fontsize', 8, 'color', 'w', ...
    'fontname', 'Times New Roman', 'FontWeight', 'bold');
m_text(111.53, 16.51, 'Ganquan Island', 'fontsize', 8, 'color', 'w', ...
    'fontname', 'Times New Roman', 'FontWeight', 'bold');
m_text(111.56, 16.55, 'Shanhu Island', 'fontsize', 8, 'color', 'w', ...
    'fontname', 'Times New Roman', 'FontWeight', 'bold');
m_text(111.642, 16.59, 'Quanfu Island', 'fontsize', 8, 'color', 'w', ...
    'fontname', 'Times New Roman', 'FontWeight', 'bold');
m_text(111.71, 16.605, 'Yin Yu', 'fontsize', 8, 'color', 'w', ...
    'fontname', 'Times New Roman', 'FontWeight', 'bold');
m_text(111.77, 16.54, 'Shi Yu', 'fontsize', 8, 'color', 'w', ...
    'fontname', 'Times New Roman', 'FontWeight', 'bold');
m_text(111.755, 16.46, 'Jinqing Island', 'fontsize', 8, 'color', 'w', ...
    'fontname', 'Times New Roman', 'FontWeight', 'bold');
m_text(111.69, 16.435, 'Chenhang Island', 'fontsize', 8, 'color', 'w', ...
    'fontname', 'Times New Roman', 'FontWeight', 'bold');
m_text(111.63, 16.435, 'Kuangzai Shoal', 'fontsize', 8, 'color', 'w', ...
    'fontname', 'Times New Roman', 'FontWeight', 'bold');

%%添加各个水道名称
m_text(111.58, 16.498, char(9312), 'fontsize', 8, 'color', 'w', 'FontWeight', 'bold');
m_text(111.59, 16.525, char(9313), 'fontsize', 8, 'color', 'w', ...
     'FontWeight', 'bold');
m_text(111.642, 16.558, char(9314), 'fontsize', 8, 'color', 'w', ...
     'FontWeight', 'bold');
m_text(111.68, 16.58, char(9315), 'fontsize', 8, 'color', 'w', ...
     'FontWeight', 'bold');
m_text(111.73, 16.56, char(9316), 'fontsize', 8, 'color', 'w', ...
   'FontWeight', 'bold');
m_text(111.725, 16.455, char(9317), 'fontsize', 8, 'color', 'w', ...
    'FontWeight', 'bold');
%%%%%标注  lagon
m_text(111.65, 16.5, 'Lagon', 'fontsize', 16, 'color', 'w', 'FontWeight', 'bold');



%% 子图 B：反射率与对数比
set(subplot(1,2,2), 'Position',[0.53, 0.10, 0.41, 0.85]);


hold on;
% 左 Y 轴：反射率
yyaxis left;
ylim([0 0.55]);
yl = ylim;
patch([111.57,111.603,111.603,111.57], [yl(1),yl(1),0.45,0.45], ...
      [0.6 0.6 0.6], 'FaceAlpha',0.3,'EdgeColor','none');
patch([111.736,111.763,111.763,111.736], [yl(1),yl(1),0.45,0.45], ...
      [0.6 0.6 0.6], 'FaceAlpha',0.3,'EdgeColor','none');
h1 = plot(lon, Band_blue(1750,:)./10000,  '-b','LineWidth',1.5);
h2 = plot(lon, Band_green(1750,:)./10000,'-g','LineWidth',1.5);
ylabel('Reflectance','FontName','Times New Roman','FontSize',14);
ax = gca; ax.YAxis(1).Color = [0 0 1];

% 右 Y 轴：对数比
yyaxis right;
xlim([min(longitude) max(longitude)]);
ylim([0.5 1.3]);
h3 = plot(lon, ln_b(1750,:)./ln_g(1750,:), '-k','LineWidth',1.5);
ylabel('Log Ratio (ln B / ln G)', ...
       'FontName','Times New Roman','FontSize',14);
ax.YAxis(2).Color = [0 0 0];

% X 轴格式
xticks2 = 111.5:0.1:111.8;
set(ax,'XTick',xticks2, ...
       'XTickLabel',arrayfun(@(x) sprintf('%.2f°E',x), xticks2,'UniformOutput',false));

% 图例和网格
legend([h1,h2,h3], {'Blue Reflectance','Green Reflectance','Log Ratio'}, ...
       'Location','northeast','FontName','Times New Roman','FontSize',16);
grid on;
ax.Box      = 'on';
ax.FontName= 'Times New Roman';
ax.FontSize= 14;
text(111.48, 0.75, 'Deep Water', 'fontsize', 16, 'color', 'k', ...
    'fontname', 'Times New Roman', 'FontWeight', 'bold');hold on
text(111.66, 0.75, 'Lagon', 'fontsize', 16, 'color', 'k', ...
    'fontname', 'Times New Roman', 'FontWeight', 'bold');       
text('Units','normalized','Position',[0.02,0.95],  'String','(B)'...
    , 'FontSize', 20, 'Color', 'k', ...
     'FontWeight', 'bold', 'BackgroundColor', 'w', ...
     'Margin', 1, 'EdgeColor', 'none');








X_bg=(ln_b./ln_g);
X_br=(ln_b./ln_r);
X_bn=(ln_b./ln_nir);
%由于激光数据存在多Nan值，提取有值区域做样本集
not_nan_T=find(~isnan(T));
T_use_select=T(not_nan_T);
% Lat_resample=Lat_resample(not_nan_T);
X_bg=X_bg(not_nan_T);  X_br=X_br(not_nan_T);  X_bn=X_bn(not_nan_T); 
X_b=Band_blue(not_nan_T);X_g=Band_green(not_nan_T);X_r=Band_red(not_nan_T); X_nir=Band_nir(not_nan_T);
not_nan_X=find(~isnan(X_bg));
% Lat_resample=Lat_resample(not_nan_X);
X_bg=X_bg(not_nan_X);  X_br=X_br(not_nan_X);  X_bn=X_bn(not_nan_X);
X_b=X_b(not_nan_X);X_g=X_g(not_nan_X);X_r=X_r(not_nan_X); X_nir=X_nir(not_nan_X);
T_use=T_use_select(not_nan_X);
data1 = [T_use,X_bg,X_br,X_bn]; 
data=data1((data1(:,1)<0) &(data1(:,1)>-20) ,:);
corrcoef(data(:,1),data(:,2))
max(data(:,2))
% 分离Y和X
Y = data(:, 1);
X = data(:, 2:4);
XX= data(:, 2);

n = length(XX);
indices = randperm(n);
% 划分训练集和测试集
train_ratio = 1;
n_train = round(train_ratio * n);
train_indices = indices(1:n_train);
test_indices = indices(n_train+1:end);
% 训练集
X_train = X(train_indices, :);
XX_train = XX(train_indices, :);
Y_train = Y(train_indices);



%% 双波段比值二次多项式
    tt=polyfit(XX_train,Y_train,2);
  waterinv1=tt(1)*XX_train.^2+tt(2)*XX_train+tt(3);

% 计算平均绝对误差（MAE）
MAE_train = mean(abs(Y_train - waterinv1));
% 计算均方根误差（RMSE）
RMSE_train = sqrt(mean((Y_train - waterinv1).^2));
% 计算拟合R^2值
R2_train = 1 - sum((Y_train - waterinv1).^2) / sum((Y_train - mean(Y_train)).^2);

h=figure(11)
x=0:0.1:20;
y1=x;
plot(x,y1,'-',-Y_train,-waterinv1,'.k');
xlim([0 20]);
ylim([0 20]);
xticks([0:2:20]);
yticks([0:2:20]);
xlabel('ICESat-2 bathymetric depth','Fontsize',14);
ylabel('Estimated depth','Fontsize',14);
txt1 = (['MAE=',num2str(MAE_train,'%.2f'),'m'])
text(14,5,txt1,'FontSize',10)
txt2 = (['RMSE=',num2str(RMSE_train,'%.2f'),'m'])
text(14,4,txt2,'FontSize',10)
txt3 = (['R^2=',num2str(R2_train,'%.2f')])
text(14,3,txt3,'FontSize',10)
txt4 = (['y=x'])
text(20,20,txt4,'FontSize',12)
text(1,19,dataname,'FontSize',10)
title('Dual-band radio quadratic polynomial model ')
picturename=('Dual-band radio quadratic polynomial model ')

cd L:\pictures\Inversion_results\20190224\Inv_Icesat2
% print(h,'-dpng','-r600',picturename);
%% 多波段二次多项式
x1=X_train(:,1);
x2=X_train(:,2);
x3=X_train(:,3);
x=[ones(length(x1),1) x1 x2 x3 (x1.^2) (x2.^2) (x3.^2)];
[b,bint,r,rint,stats]=regress(Y_train,x)
waterinv2=x*b;

% 计算平均绝对误差（MAE）
MAE_train = mean(abs(Y_train - waterinv2));
% 计算均方根误差（RMSE）
RMSE_train = sqrt(mean((Y_train - waterinv2).^2));
% 计算拟合R^2值
R2_train = 1 - sum((Y_train - waterinv2).^2) / sum((Y_train - mean(Y_train)).^2);
% 显示结果

fprintf('训练集平均绝对误差（MAE）：%.4f\n', MAE_train);
fprintf('训练集均方根误差（RMSE）：%.4f\n', RMSE_train);
fprintf('训练集拟合R^2值：%.4f\n', R2_train);
h= figure(12)
x=0:0.1:20;
y1=x;
plot(x,y1,'-',-Y_train,-waterinv1,'.k');
xlim([0 20]);
ylim([0 20]);
xticks([0:2:20])
yticks([0:2:20])
xlabel('ICESat-2 bathymetric depth','Fontsize',14);
ylabel('Estimated depth','Fontsize',14);
txt1 = (['MAE=',num2str(MAE_train,'%.2f'),'m'])
text(14,5,txt1,'FontSize',10)
txt2 = (['RMSE=',num2str(RMSE_train,'%.2f'),'m'])
text(14,4,txt2,'FontSize',10)
txt3 = (['R^2=',num2str(R2_train,'%.2f')])
text(14,3,txt3,'FontSize',10)
txt4 = (['y=x'])
text(20,20,txt4,'FontSize',12)
text(1,19,dataname,'FontSize',10)
title('Multiband radio quadratic polynomial model ')
picturename=('Multiband radio quadratic polynomial model ')
% print(h,'-dpng','-r600',picturename);

%% y=aX1+b 双波段对数比值模型
XX_train = XX(train_indices, :);
XX_train = [XX_train, ones(size(XX_train, 1), 1)];
% 拟合线性回归模型
beta_2bands = (XX_train' * XX_train) \ (XX_train' * Y_train);
% 显示拟合方程
fprintf('拟合方程: Y = %.4f*X1 + %.4f\n', beta_2bands(1), beta_2bands(2));
% 训练集预测
Y_train_pred = XX_train * beta_2bands;
% 计算平均相对误差（MRE）
MRE_train = mean(abs((Y_train - Y_train_pred) ./ Y_train));
% 计算平均绝对误差（MAE）
MAE_train = mean(abs(Y_train - Y_train_pred));
% 计算均方根误差（RMSE）
RMSE_train = sqrt(mean((Y_train - Y_train_pred).^2));
% 计算拟合R^2值
R2_train = 1 - sum((Y_train - Y_train_pred).^2) / sum((Y_train - mean(Y_train)).^2);
% 显示结果
fprintf('训练集平均相对误差（MRE）：%.4f\n', MRE_train);
fprintf('训练集平均绝对误差（MAE）：%.4f\n', MAE_train);
fprintf('训练集均方根误差（RMSE）：%.4f\n', RMSE_train);
fprintf('训练集拟合R^2值：%.4f\n', R2_train);

h=figure(13)
x=0:0.1:20;
y1=x;
plot(x,y1,'-',-Y_train,-Y_train_pred,'.k');
xlim([0 20]);
ylim([0 20]);
xticks([0:2:20])
yticks([0:2:20])
xlabel('ICESat-2 bathymetric depth','Fontsize',14);
ylabel('Estimated depth','Fontsize',14);
txt1 = (['MAE=',num2str(MAE_train,'%.2f'),'m'])
text(14,5,txt1,'FontSize',10)
txt2 = (['RMSE=',num2str(RMSE_train,'%.2f'),'m'])
text(14,4,txt2,'FontSize',10)
txt3 = (['R^2=',num2str(R2_train,'%.2f')])
text(14,3,txt3,'FontSize',10)
txt4 = (['y=x'])
text(20,20,txt4,'FontSize',12)
text(1,19,dataname,'FontSize',10)
% txt5 = (['Y = ',num2str(beta_2bands(1),'%.2f'),'*X1 +' ,num2str(beta_2bands(2),'%.2f')])
% text(14,2,txt5,'FontSize',10)
text(1,19,dataname,'FontSize',10)
title('Dual-band radio model ')

picturename=('Dual-band radio model ')
% print(h,'-dpng','-r600',picturename);

%% y=aX1+bX2+cX3+d多波段对数比值模型
% 添加常数项
x = [X_train, ones(size(X_train, 1), 1)];
Y=Y_train;
% 拟合线性回归模型
beta_multi_band = (x' * x) \ (x' * Y);
% 显示拟合方程
fprintf('拟合方程: Y = %.4f*X1 + %.4f*X2 + %.4f*X3 + %.4f\n', beta_multi_band(1), beta_multi_band(2), beta_multi_band(3), beta_multi_band(4));
% 训练集预测
Y_train_pred = x * beta_multi_band;
% 计算平均相对误差（MRE）
MRE_train = mean(abs((Y_train - Y_train_pred) ./ Y_train));
% 计算平均绝对误差（MAE）
MAE_train = mean(abs(Y_train - Y_train_pred));
% 计算均方根误差（RMSE）
RMSE_train = sqrt(mean((Y_train - Y_train_pred).^2));
% 计算拟合R^2值
R2_train = 1 - sum((Y_train - Y_train_pred).^2) / sum((Y_train - mean(Y_train)).^2);
% 显示结果
fprintf('训练集平均相对误差（MRE）：%.4f\n', MRE_train);
fprintf('训练集平均绝对误差（MAE）：%.4f\n', MAE_train);
fprintf('训练集均方根误差（RMSE）：%.4f\n', RMSE_train);
fprintf('训练集拟合R^2值：%.4f\n', R2_train);

h=figure(14)
% set(gcf,'position',[left,bottom,width,height]);
% set(gcf,'position',[400 150 600 500])
x=0:0.1:20;
y1=x;
plot(x,y1,'-',-Y_train,-Y_train_pred,'.k','linewid',1, 'MarkerSize', 7);
box on
xticks([0:4:20])
yticks([0:4:20])
ax = gca;
set( ax, 'xlim', [0 20], 'fontsize', 18)
set( ax, 'ylim', [0 20], 'fontsize', 18)
xlabel('ICESat-2 bathymetric depth(m)','Fontsize',18);
ylabel('Estimated depth(m)','Fontsize',18);
% title('Multiband radio model ')
txt1 = (['MAE=',num2str(MAE_train,'%.2f'),'m'])
text(12,5.5,txt1,'FontSize',18)
txt2 = (['RMSE=',num2str(RMSE_train,'%.2f'),'m'])
text(12,3.5,txt2,'FontSize',18)
R=sqrt(R2_train);
txt3 = (['R=',num2str(R,'%.2f')])
text(12,1.5,txt3,'FontSize',18)
txt4 = (['y=x'])
text(20,20,txt4,'FontSize',18,'color','b')
% text(1,19,dataname,'FontSize',14)

picturename=('Multiband radio model ')
% print(h,'-dpng','-r600',picturename);
%% 区域水深反演
%利用反演公式结合反射率数据集 反演珊瑚岛区域水深值
fanshelv(:,7)=log(fanshelv(:,3))./log(fanshelv(:,4));%%反射率因子放大10000倍
fanshelv(:,8)=log(fanshelv(:,3))./log(fanshelv(:,5));%%反射率因子放大10000倍
fanshelv(:,9)=log(fanshelv(:,3))./log(fanshelv(:,6));%%反射率因子放大10000倍

X1=fanshelv(:,7);
X2=fanshelv(:,8);
X3=fanshelv(:,9);
xx=[(X1.^2) X1 ones(length(X1),1)];
xxx=[ones(length(X1),1) X1 X2 X3 (X1.^2) (X2.^2) (X3.^2)];

SDB_2Band_radio=(beta_2bands(1)).*X1+beta_2bands(2);
SDB_multi_Band_radio=(beta_multi_band(1)).*X1+beta_multi_band(2)*X2+beta_multi_band(3)*X3+beta_multi_band(4);
SDB_2Band_radio_quadratic_polynomial=xx*tt';
SDB_multi_Band_radio_quadratic_polynomial=xxx*b;


fanshelv(:,10)=SDB_2Band_radio;
fanshelv(:,11)=SDB_multi_Band_radio;
fanshelv(:,12)=SDB_2Band_radio_quadratic_polynomial;
fanshelv(:,13)=SDB_multi_Band_radio_quadratic_polynomial;
%% 绘制包含陆地和伪浅海的永乐环礁地形图
% situ_Depth1=situ_Depth';
% % imshow(situ_Depth1,[])
% situ_Depth2=situ_Depth1(:);

% % %%%%多波段比值法
% figure(444)
% % set(gcf,'position',[left,bottom,width,height]);
% set(gcf,'position',[400 100 1100 600])
% set(gca,'position', [0.10 0.1 0.79 0.86]);
% % 设置投影和边界框
% m_proj('Mercator', 'longitudes', [min(fanshelv(:, 2)) max(fanshelv(:, 2))], 'latitudes', [min(fanshelv(:, 1)) max(fanshelv(:, 1))]);%%%甘泉岛范围
% % 绘制散点图
% m_scatter(fanshelv(:, 2), fanshelv(:, 1), 2, fanshelv(:,13), 'filled'); % 15是点的大小，可以根据需要调整
% custom_colormap = jet(10);  % 可以选择其他颜色图，如 parula(10), hot(10) 等
% % 应用颜色映射
% colormap(custom_colormap);
% caxis([-20 0]);
% colorbar('ylim',[-20,0],'ytick',[-20,-16,-12,-8,-4,0]); % 添加颜色条，表示水深
% h = colorbar;
% h.Label.String = 'Depth(m)';
% set(h,'FontSize',16);
% % 添加经纬度边框和网格
% m_grid('box', 'on', 'tickdir', 'in','linewi',1.5,'fontsize', 20, 'tickstyle', 'dd');
% % % 添加指北针
% % h = m_northarrow(111.485, 16.57, 0.02, 'type', 2, 'linewi', 1, 'facecolor', 'K', 'aspect', 1.09);
% % % %添加比例尺
% % m_ruler([0.13 0.33],0.73,4,'color','b',...
% %     'linewid',5,'tickdir','out','ticklen',0.01,'fontsize',20);
% % 设置经纬度标签
% % set(gca, 'XTickLabel', num2str(get(gca, 'XTick')'));
% % set(gca, 'YTickLabel', num2str(get(gca, 'YTick')'));
% % txt1 = (['Date: 2022-03-15'])
% % m_text(111.48,16.61,txt1,'FontSize',22)
% xlabel('Longitude (deg)','Fontsize',20);
% ylabel('Latitude (deg)','Fontsize',20);

%% 水陆分离
% %计算NDWI
green_band=fanshelv(:,4);
nir_band=fanshelv(:,6);
ndwi = (green_band - nir_band) ./ (green_band + nir_band);
%设定阈值进行水陆分离（一般选择0作为阈值）
threshold = 0;
Land_mask = ndwi < threshold;
% 提取数据
latitude = fanshelv(:, 1); % 第一列为纬度
longitude = fanshelv(:, 2); % 第二列为经度
Land_mask_longitude=longitude(Land_mask);
Land_mask_latitude=latitude(Land_mask);
Land_mask=[Land_mask_longitude,Land_mask_latitude];

depth_2bands = fanshelv(:, 10); % 第10列为水深
depth_2bands(fanshelv(:,3)<690)=nan;
depth_2bands (Land_mask==1) = nan;%%剔除反演大于0的值，视为陆地

depth_multi_band = fanshelv(:, 11); 
depth_multi_band(fanshelv(:,3)<690)=nan;
depth_multi_band (Land_mask==1) = nan;%%剔除反演大于0的值，视为陆地

SDB_2Band_radio_quadratic_polynomial(fanshelv(:,3)>690)=nan;
SDB_2Band_radio_quadratic_polynomial (Land_mask==1) = nan;%%剔除反演大于0的值，视为陆地

SDB_multi_Band_radio_quadratic_polynomial(fanshelv(:,3)<690)=nan;
SDB_multi_Band_radio_quadratic_polynomial (Land_mask==1) = nan;%%剔除反演大于0的值，视为陆地
%% % %% 绘制永乐环礁地形图
situ_Depth1=situ_Depth';
% imshow(situ_Depth1,[])
situ_Depth2=situ_Depth1(:);

% %%%%多波段比值法
figure(333)
% set(gcf,'position',[left,bottom,width,height]);
set(gcf,'position',[400 100 1100 600])
set(gca,'position', [0.10 0.1 0.79 0.86]);
% 设置投影和边界框
m_proj('Mercator', 'longitudes', [min(longitude) max(longitude)], 'latitudes', [min(latitude) max(latitude)]);%%%甘泉岛范围
% 绘制散点图
m_scatter(longitude, latitude, 2, depth_multi_band, 'filled'); % 15是点的大小，可以根据需要调整
custom_colormap = jet(10);  % 可以选择其他颜色图，如 parula(10), hot(10) 等
% 应用颜色映射
colormap(custom_colormap);
caxis([-20 0]);
hold on
%%绘制陆地（提取Z为NaN的点，绘制为黑色）
load('L:\code\Sentienl2_Bathymetry\Land_mask.mat');
m_scatter(Land_mask(:,1), Land_mask(:,2), 1, [0.2 0.2 0.2], 'filled'); % 黑色

colorbar('ylim',[-20,0],'ytick',[-20,-16,-12,-8,-4,0]); % 添加颜色条，表示水深
h = colorbar;
h.Label.String = 'Depth(m)';
set(h,'FontSize',16);
% 添加经纬度边框和网格
m_grid('box', 'on', 'tickdir', 'in','linewi',1.5,'fontsize', 20, 'tickstyle', 'dd');
% 添加指北针
h = m_northarrow(111.485, 16.57, 0.02, 'type', 2, 'linewi', 1, 'facecolor', 'K', 'aspect', 1.09);
% %添加比例尺
m_ruler([0.13 0.33],0.73,4,'color','b',...
    'linewid',5,'tickdir','out','ticklen',0.01,'fontsize',20);

txt1 = (['Date: 2022-03-15'])
m_text(111.48,16.61,txt1,'FontSize',22)
xlabel('Longitude (deg)','Fontsize',20);
ylabel('Latitude (deg)','Fontsize',20);

picturename=('Yongle Atoll Bathymetry_20220315')

% print(gcf,'-dpng','-r600',picturename);
% % 保存文件
% filename4 = [ 'Yongle_Atoll_SDB_coverland.txt'];
% filename4 = [ 'Yongle_Atoll_SDB.txt'];
% 
% filename_full = ['L:\code\Sentienl2_Bathymetry\result\', filename4];
% 
% SDB_result=[longitude,latitude,depth_multi_band];
% % 打开文件以写入模式
% fileID = fopen(filename_full,'w');
% if fileID == -1
%     error('文件无法打开');
% end
% formatSpec = '%.16f,%.16f,%.16f\n';
% % 写入表头
% fprintf(fileID, '经度,纬度,水深\n');
% % 写入数据
% for j = 1:size(SDB_result, 1)
%     fprintf(fileID, formatSpec, double(SDB_result(j, 1)), double(SDB_result(j, 2)), double(-SDB_result(j, 3)));
% end
% % 关闭文件
% fclose(fileID);
% disp(['数据已保存到 ', filename_full]);

%% 
%验证数据存在多Nan值，提取有值区域做样本集
depth_2bands=(reshape(depth_2bands,cols,rows))';
depth_multi_band=(reshape(depth_multi_band,cols,rows))';

SDB_2Band_radio_quadratic_polynomial=(reshape(SDB_2Band_radio_quadratic_polynomial,cols,rows))';
SDB_multi_Band_radio_quadratic_polynomial=(reshape(SDB_multi_Band_radio_quadratic_polynomial,cols,rows))';

not_nan_situ_Depth=find(~isnan(situ_Depth));
situ_Depth_select1=situ_Depth(not_nan_situ_Depth);%%永乐环礁区域重采样的原位水深点
depth_2band_select1=depth_2bands(not_nan_situ_Depth);%%永乐环礁区域重采样的原位水深点处的反演水深值（双波段对数比值）
depth_multi_band_select1=depth_multi_band(not_nan_situ_Depth);%%永乐环礁区域重采样的原位水深点处的反演水深值（多波段对数比值法）
SDB_2Band_radio_quadratic_polynomial_select1=SDB_2Band_radio_quadratic_polynomial(not_nan_situ_Depth);
SDB_multi_Band_radio_quadratic_polynomial_select1=SDB_multi_Band_radio_quadratic_polynomial(not_nan_situ_Depth);
%由于将陆地和截止反射率位置的水深值置为了nan,需要筛掉
not_nan_depth_2bands_select1=find(~isnan(depth_2band_select1));
% not_nan_depth_multi_band_select1=find(~isnan(depth_multi_band_select1));%%两种反演方法非nan值区域是一样的
depth_2bands_select2=depth_2band_select1(not_nan_depth_2bands_select1);
depth_multi_band_select2=depth_multi_band_select1(not_nan_depth_2bands_select1);
situ_Depth_select2=situ_Depth_select1(not_nan_depth_2bands_select1);%%%永乐环礁区域重采样的原位水深点
SDB_2Band_radio_quadratic_polynomial_select2=SDB_2Band_radio_quadratic_polynomial_select1(not_nan_depth_2bands_select1);
SDB_multi_Band_radio_quadratic_polynomial_select2=SDB_multi_Band_radio_quadratic_polynomial_select1(not_nan_depth_2bands_select1);
%%%%%筛选出原位水深数据中0-20米的数据
mmm=find(situ_Depth_select2<20);
situ_Depth_select3=situ_Depth_select2(mmm);
depth_2bands_select3=depth_2bands_select2(mmm);
depth_multi_band_select3=depth_multi_band_select2(mmm);
SDB_2Band_radio_quadratic_polynomial_select3=SDB_2Band_radio_quadratic_polynomial_select2(mmm);
SDB_multi_Band_radio_quadratic_polynomial_select3=SDB_multi_Band_radio_quadratic_polynomial_select2(mmm);


%%%%%绘图  Sentienl2反演水深与原位水深
%% 双波段比值反演
% 计算平均绝对误差（MAE）
MAE_train = mean(abs(situ_Depth_select3 - depth_2bands_select3));
% 计算均方根误差（RMSE）
RMSE_train = sqrt(mean((situ_Depth_select3 - depth_2bands_select3).^2));
% 计算拟合R^2值
R2_train = 1 - sum((situ_Depth_select3 - depth_2bands_select3).^2) / sum((situ_Depth_select3 - mean(situ_Depth_select3)).^2);
% 显示结果
fprintf('训练集平均绝对误差（MAE）：%.4f\n', MAE_train);
fprintf('训练集均方根误差（RMSE）：%.4f\n', RMSE_train);
fprintf('训练集拟合R^2值：%.4f\n', R2_train);

depth_2bands_select3=-depth_2bands_select3;
% 计算平均绝对误差（MAE）
MAE_train = mean(abs(situ_Depth_select3 - depth_2bands_select3));
% 计算均方根误差（RMSE）
RMSE_train = sqrt(mean((situ_Depth_select3 - depth_2bands_select3).^2));
% 计算拟合R^2值
R2_train = 1 - sum((situ_Depth_select3 - depth_2bands_select3).^2) / sum((situ_Depth_select3 - mean(situ_Depth_select3)).^2);
% 显示结果
fprintf('训练集平均绝对误差（MAE）：%.4f\n', MAE_train);
fprintf('训练集均方根误差（RMSE）：%.4f\n', RMSE_train);
fprintf('训练集拟合R^2值：%.4f\n', R2_train);
h=figure(21)
x=0:0.1:20;
y1=x;
plot(x,y1,'-',situ_Depth_select3,depth_2bands_select3,'.k');
xlim([0 20]);
ylim([0 20]);
xticks([0:2:20])
yticks([0:2:20])
xlabel('Si-tu depth depth','Fontsize',14);
ylabel('Estimated depth','Fontsize',14);
txt1 = (['MAE=',num2str(MAE_train,'%.2f'),'m'])
text(14,5,txt1,'FontSize',10)
txt2 = (['RMSE=',num2str(RMSE_train,'%.2f'),'m'])
text(14,4,txt2,'FontSize',10)
txt3 = (['R^2=',num2str(R2_train,'%.2f')])
text(14,3,txt3,'FontSize',10)
txt4 = (['y=x'])
text(20,20,txt4,'FontSize',12)
text(1,19,dataname,'FontSize',10)
title('Dual-band radio model ')


cd L:\pictures\Inversion_results\20190224/Inv_situ
picturename=('Dual-band radio model ')
% print(h,'-dpng','-r600',picturename);
%% 多波段比值反演
% 计算平均绝对误差（MAE）
MAE_train = mean(abs(situ_Depth_select3 - depth_multi_band_select3));
% 计算均方根误差（RMSE）
RMSE_train = sqrt(mean((situ_Depth_select3 - depth_multi_band_select3).^2));
% 计算拟合R^2值
R2_train = 1 - sum((situ_Depth_select3 - depth_multi_band_select3).^2) / sum((situ_Depth_select3 - mean(situ_Depth_select3)).^2);
% 显示结果
fprintf('训练集平均绝对误差（MAE）：%.4f\n', MAE_train);
fprintf('训练集均方根误差（RMSE）：%.4f\n', RMSE_train);
fprintf('训练集拟合R^2值：%.4f\n', R2_train);

depth_multi_band_select3=-depth_multi_band_select3;
% 计算平均绝对误差（MAE）
MAE_train = mean(abs(situ_Depth_select3 - depth_multi_band_select3));
% 计算均方根误差（RMSE）
RMSE_train = sqrt(mean((situ_Depth_select3 - depth_multi_band_select3).^2));
% 计算拟合R^2值
R2_train = 1 - sum((situ_Depth_select3 - depth_multi_band_select3).^2) / sum((situ_Depth_select3 - mean(situ_Depth_select3)).^2);
% 显示结果
fprintf('训练集平均绝对误差（MAE）：%.4f\n', MAE_train);
fprintf('训练集均方根误差（RMSE）：%.4f\n', RMSE_train);
fprintf('训练集拟合R^2值：%.4f\n', R2_train);



h=figure(22)
x=0:0.1:20;
y1=x;
plot(x,y1,'-',situ_Depth_select3,depth_multi_band_select3,'.k','linewid',1, 'MarkerSize', 7);
xticks([0:4:20])
yticks([0:4:20])
ax = gca;
set( ax, 'xlim', [0 20], 'fontsize', 18)
set( ax, 'ylim', [0 20], 'fontsize', 18)
xlabel('Si-tu depth(m)','Fontsize',18);
ylabel('Estimated depth(m)','Fontsize',18);
% xlabel('ICESat-2 bathymetric depth(m)','Fontsize',18);
% ylabel('Estimated depth(m)','Fontsize',18);
% title('Multiband radio model ')
txt1 = (['MAE=',num2str(MAE_train,'%.2f'),'m'])
text(12,5.5,txt1,'FontSize',18)
txt2 = (['RMSE=',num2str(RMSE_train,'%.2f'),'m'])
text(12,3.5,txt2,'FontSize',18)
txt3 = (['R^2=',num2str(R2_train,'%.2f')])
text(12,1.5,txt3,'FontSize',18)
txt4 = (['y=x'])
text(20,20,txt4,'FontSize',18,'color','b')
box on
% print(h,'-dpng','-r600',picturename);








%% 双波段比值二次项反演
% 计算平均绝对误差（MAE）
MAE_train = mean(abs(situ_Depth_select3 - SDB_2Band_radio_quadratic_polynomial_select3));
% 计算均方根误差（RMSE）
RMSE_train = sqrt(mean((situ_Depth_select3 - SDB_2Band_radio_quadratic_polynomial_select3).^2));
% 计算拟合R^2值
R2_train = 1 - sum((situ_Depth_select3 - SDB_2Band_radio_quadratic_polynomial_select3).^2) / sum((situ_Depth_select3 - mean(situ_Depth_select3)).^2);
% 显示结果
fprintf('训练集平均绝对误差（MAE）：%.4f\n', MAE_train);
fprintf('训练集均方根误差（RMSE）：%.4f\n', RMSE_train);
fprintf('训练集拟合R^2值：%.4f\n', R2_train);

SDB_2Band_radio_quadratic_polynomial_select3=-SDB_2Band_radio_quadratic_polynomial_select3;
% 计算平均绝对误差（MAE）
MAE_train = mean(abs(situ_Depth_select3 - SDB_2Band_radio_quadratic_polynomial_select3));
% 计算均方根误差（RMSE）
RMSE_train = sqrt(mean((situ_Depth_select3 - SDB_2Band_radio_quadratic_polynomial_select3).^2));
% 计算拟合R^2值
R2_train = 1 - sum((situ_Depth_select3 - SDB_2Band_radio_quadratic_polynomial_select3).^2) / sum((situ_Depth_select3 - mean(situ_Depth_select3)).^2);
% 显示结果
fprintf('训练集平均绝对误差（MAE）：%.4f\n', MAE_train);
fprintf('训练集均方根误差（RMSE）：%.4f\n', RMSE_train);
fprintf('训练集拟合R^2值：%.4f\n', R2_train);

h=figure(23)
x=0:0.1:20;
y1=x;
plot(x,y1,'-',situ_Depth_select3,SDB_2Band_radio_quadratic_polynomial_select3,'.k');
xlim([0 20]);
ylim([0 20]);
xticks([0:2:20])
yticks([0:2:20])
xlabel('Si-tu depth','Fontsize',14);
ylabel('Estimated depth','Fontsize',14);
txt1 = (['MAE=',num2str(MAE_train,'%.2f'),'m'])
text(14,5,txt1,'FontSize',10)
txt2 = (['RMSE=',num2str(RMSE_train,'%.2f'),'m'])
text(14,4,txt2,'FontSize',10)
txt3 = (['R^2=',num2str(R2_train,'%.2f')])
text(14,3,txt3,'FontSize',10)
txt4 = (['y=x'])
text(20,20,txt4,'FontSize',12)
text(1,19,dataname,'FontSize',10)
title('Dual-band quadratic polynomial model')
picturename=('Dual-band quadratic polynomial model')
% print(h,'-dpng','-r600',picturename);
%% 双波段比值二次项反演
% 计算平均绝对误差（MAE）
MAE_train = mean(abs(situ_Depth_select3 - SDB_multi_Band_radio_quadratic_polynomial_select3));
% 计算均方根误差（RMSE）
RMSE_train = sqrt(mean((situ_Depth_select3 - SDB_multi_Band_radio_quadratic_polynomial_select3).^2));
% 计算拟合R^2值
R2_train = 1 - sum((situ_Depth_select3 - SDB_multi_Band_radio_quadratic_polynomial_select3).^2) / sum((situ_Depth_select3 - mean(situ_Depth_select3)).^2);
% 显示结果
fprintf('训练集平均绝对误差（MAE）：%.4f\n', MAE_train);
fprintf('训练集均方根误差（RMSE）：%.4f\n', RMSE_train);
fprintf('训练集拟合R^2值：%.4f\n', R2_train);

SDB_multi_Band_radio_quadratic_polynomial_select3=-SDB_multi_Band_radio_quadratic_polynomial_select3;
% 计算平均绝对误差（MAE）
MAE_train = mean(abs(situ_Depth_select3 - SDB_multi_Band_radio_quadratic_polynomial_select3));
% 计算均方根误差（RMSE）
RMSE_train = sqrt(mean((situ_Depth_select3 - SDB_multi_Band_radio_quadratic_polynomial_select3).^2));
% 计算拟合R^2值
R2_train = 1 - sum((situ_Depth_select3 - SDB_multi_Band_radio_quadratic_polynomial_select3).^2) / sum((situ_Depth_select3 - mean(situ_Depth_select3)).^2);
% 显示结果
fprintf('训练集平均绝对误差（MAE）：%.4f\n', MAE_train);
fprintf('训练集均方根误差（RMSE）：%.4f\n', RMSE_train);
fprintf('训练集拟合R^2值：%.4f\n', R2_train);

h=figure(24)
x=0:0.1:20;
y1=x;
plot(x,y1,'-',situ_Depth_select3,SDB_multi_Band_radio_quadratic_polynomial_select3,'.k','linewid',1, 'MarkerSize', 7);
xlim([0 20]);
xticks([0:4:20])
yticks([0:4:20])
ax = gca;
set( ax, 'xlim', [0 20], 'fontsize', 18)
set( ax, 'ylim', [0 20], 'fontsize', 18)
xlabel('Measured depth(m)','Fontsize',18);
ylabel('Estimated depth(m)','Fontsize',18);
% xlabel('ICESat-2 bathymetric depth(m)','Fontsize',18);
% ylabel('Estimated depth(m)','Fontsize',18);
% title('Multiband radio model ')
txt1 = (['MAE=',num2str(MAE_train,'%.2f'),'m'])
text(12,5.5,txt1,'FontSize',18)
txt2 = (['RMSE=',num2str(RMSE_train,'%.2f'),'m'])
text(12,3.5,txt2,'FontSize',18)
R=sqrt(R2_train);
txt3 = (['R=',num2str(R,'%.2f')])
text(12,1.5,txt3,'FontSize',18)
txt4 = (['y=x'])
text(20,20,txt4,'FontSize',18,'color','b')
box on
picturename=('Multiband quadratic polynomial model')

% print(h,'-dpng','-r600',picturename);


h=figure(22)
x=0:0.1:20;
y1=x;
plot(x,y1,'-',situ_Depth_select3,-depth_multi_band_select3,'.k','linewid',1, 'MarkerSize', 7);
xticks([0:4:20])
yticks([0:4:20])
ax = gca;
set( ax, 'xlim', [0 20], 'fontsize', 18)
set( ax, 'ylim', [0 20], 'fontsize', 18)
xlabel('Si-tu depth(m)','Fontsize',18);
ylabel('Estimated depth(m)','Fontsize',18);
% xlabel('ICESat-2 bathymetric depth(m)','Fontsize',18);
% ylabel('Estimated depth(m)','Fontsize',18);
% title('Multiband radio model ')
txt1 = (['MAE=',num2str(MAE_train,'%.2f'),'m'])
text(12,5.5,txt1,'FontSize',18)
txt2 = (['RMSE=',num2str(RMSE_train,'%.2f'),'m'])
text(12,3.5,txt2,'FontSize',18)
txt3 = (['R^2=',num2str(R2_train,'%.2f')])
text(12,1.5,txt3,'FontSize',18)
txt4 = (['y=x'])
text(20,20,txt4,'FontSize',18,'color','b')
box on
