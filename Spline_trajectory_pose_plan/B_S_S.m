% 清空环境
clear; clc;
close all;
% 定义控制点
% P = [1 1 1; 2 3 2; 4 5 5; 2 3 3; 5 4 3; 6 7 1; 9 9 8; 12 15 11]; % 控制点坐标
% P = [1 1 1; 2 3 2; 4 5 5; 6 6 3; 5 5 1; 6 7 1; 9 9 8; 12 15 11];
% P=[1623.382813, 1287.131104, 1209.996582, 1065.779297, 983.905090, 892.094543, 835.077087, 796.471741, 751.769348,713.332092,669.832581,601.487244,556.960999,475.868347,418.927307,344.341064, 300.111877, 253.929565, 200.157654, 161.352356, 106.556396, 49.665283, -6.967102
%    704.965027, 689.627258, 694.901245, 700.669373, 704.640991,707.915222, 709.012573, 703.724487, 698.676758,692.176392,690.623779,690.779785,690.075928,691.352539,695.278809,699.948914, 702.718506, 688.662476,677.834961, 669.927490, 661.460449, 654.651489, 651.501038
%    -12.510620, 251.928284, 217.869019, 193.214111, 111.798584,11.280273, 8.179810, 123.177490, 215.058594,273.975037,304.385864,324.083374,296.066772,270.943787,198.506470,117.120605, 41.464844, 157.157593, 251.385193, 294.979492, 327.523560, 331.231323, 307.992493 ];
% 	
% P= P'/1000;
% P = [0,0,0,90.001/180*pi, -0.001/180*pi, 90.000/180*pi;
%     0.140765,-0.502160, 0.123521, 0.001/180*pi, -0.001/180*pi, 90.000/180*pi;
%     -0.218351,-0.474313, 0.135510, -60.391/180*pi, -0.920/180*pi, 88.689/180*pi;
%     -0.220025,-0.492161, 0.199772, -60.401/180*pi,  0.149/180*pi, 90.215/180*pi;
%      0.006761,-0.540792, 0.296533, 0.528/180*pi, -0.962/180*pi, 87.596/180*pi];
% P = [
% -0.004033,0.3055,0.3909,1.1784,1.5462,-0.83048;
% -0.002464,0.30846,0.28987,1.2859,1.5521,-0.72293;
% -0.084906,0.47651,0.29636,1.2525,1.5337,-0.76819;
% -0.020277,0.50181,0.4644,1.2534,1.5244,-0.76136;
% 0.030274,0.39687,0.46005,1.2622,1.5269,-0.74918;
% -0.006156,0.30654,0.59421,1.2058,1.5164,-0.80915
% ];

P = [
 -388 -426 195 90 0 -20
 -112 -426 195 90 0 -10
 153 -525 231 90 0 12
-202 -249 200 90 0 260
-302 -49 200 90 0 230
-227 226 200 90 0 180
117 326 196 90 0 110
117 326 86 90 0 110
];
P(:,1:3) = P(:,1:3)/1000;
P(:,4:6) = P(:,4:6)/pi*180;

%   pose = [0.0, 0.0, 0.0;
%   0.001/180*pi, -0.001/180*pi, 90.000/180*pi;
% -40.391/180*pi, -0.920/180*pi, 88.689/180*pi;
% -40.401/180*pi,  0.149/180*pi, 90.215/180*pi;
%  -4.528/180*pi, -0.962/180*pi, 87.596/180*pi];
% R = [];
% r0 = eul2rotm(pose(1,:),"ZYX");
% for i = 2:5
%     r = eul2rotm(pose(i,:),"ZYX");
%     delta = r0' * r;
%     R = [R;delta];
%     r0 = r;
% end
% eqradius = 0.1;
% [k, alpha] = rotationMatrixToAxisAngle(R(1:3,1:3));
% 
% for i = 0:0.01:1
%     rr = axisAngleToRotationMatrix(k, alpha*i);
%     aa = rotm2eul(rr,"ZYX")
% end
% if (alpha ~= 0 && alpha * eqradius > dist)
% 	pathlength = alpha * eqradius;
% 	scalerot = 1 / eqradius;
% 	scalelin = dist / pathlength;
% elseif (dist ~= 0)
% 	pathlength = dist;
% 	scalerot = alpha / pathlength;
% 	scalelin = 1;
% else
% 	pathlength = 0;
% 	scalerot = 1;
% 	scalelin = 1;
% end
% alphaLen = alpha;


n = size(P, 1) - 1;  % 控制点数量减 1
% P = P/1000;
% 定义 B 样条曲线的次数
p = 3;  % 3 次 B 样条曲线
steps = 500;
% 定义节点向量
% knots = [zeros(1, p), linspace(0, 1, n - p + 2), ones(1, p)];
knots = knots_uniform(n, p);
% 离散化参数 u
u = linspace(0, 1, steps);  % 100 个离散点

% 初始化位置、速度和加速度
C = zeros(length(u), 3);  % 位置
V = zeros(length(u), 3);  % 速度
A = zeros(length(u), 3);  % 加速度
J = zeros(length(u), 3);  % 加加速度
% bspline_basis(p, knots, 0.99999999)
% 初始化速度向量的模
speed_norm = zeros(size(u));
kappa = zeros(length(u), 1);  % 曲率半径

% 计算位置、速度和加速度
for i = 1:length(u)
    % 计算基函数及其导数
    if u(i) == 1
        u(i) = u(i) - 1e-10;
    end
    [N0, N1, N2, N3] = bspline_basis(p, knots, u(i));

    % 计算位置
    C(i, :) = N0 * P(:,1:3);
    
    % 计算速度
    V(i, :) = N1 * P(:,1:3);
    
    % 计算加速度
    A(i, :) = N2 * P(:,1:3);

    % 计算加速度
    J(i, :) = N3 * P(:,1:3);

    cross_prod = cross(V(i, :), A(i, :));  % 二维曲线需要补零
    kappa(i) = 1/(norm(cross_prod(3)) / norm(V(i, :))^3);
%     cross_prod = cross(C1, C2);  % 二维曲线需要补零
%     kappa(i) = 1 / (norm(cross_prod(3)) / norm(C1)^3);
    % 计算速度向量
    C_prime = N1 * P(:,1:3);
    
    % 计算速度向量的模
    speed_norm(i) = norm(C_prime);
end
kappa(steps) = 100;
% plot3(C(:,1), C(:,2), C(:,3));
% l = trapz(u, speed_norm)
L = [];
l_last = 0;
LL = [];
% 数值积分计算曲线长度
for i = 2:length(u)
    l = trapz(u(1:i), speed_norm(1:i));
    L = [L; l - l_last];
    LL = [LL; l];
    l_last = l;
end
l_sum = sum(L);
Vmax = 0.03;
V = Vmax;
threshold = 0.3;
% kappa = [22.0454076850486	22.1198252322794	22.1957589621591	22.2731258002004	22.3518418876507	22.4318225613703	22.5129823353024	22.5952348836060	22.6784930255216	22.7626687120534	22.8476730145500	22.9334161152806	23.0198073001071	23.1067549533633	23.1941665550632	23.2819486805718	23.3700070028824	23.4582462976632	23.5465704512469	23.6348824717576	23.7230845035891	23.8110778454691	23.8987629723701	23.9860395615542	24.0728065230711	24.1589620350620	24.2444035842618	24.3290280121357	24.4127315671350	24.4954099636131	24.5769584480044	24.6572718729406	24.7362447800590	24.8137714923465	24.8897462169700	24.9640631596575	25.0366166518318	25.1073012918505	25.1760121018806	25.2426447021412	25.3070955044752	25.3692619274805	25.4290426357413	25.4863378060526	25.5410494239509	25.5930816143441	25.6423410105963	25.6887371670843	25.7321830210124	25.7725954101850	25.8098956545062	25.8440102102485	25.8748714076418	25.9024182841321	25.9265975278088	25.9473645480862	25.9646846938307	25.9785346428956	25.9889039915931	25.9957970782177	25.9992350815710	25.9992584438623	25.9959296777847	25.9893366305327	25.9795962937574	25.9668592688636	25.9513150229035	25.9331981032383	25.9127955213775	25.8904555709564	25.8665984158333	25.8417288774727	25.8164519741088	25.7914919288944	25.7677155863621	25.7461614792020	25.7280762044258	25.7149603495890	25.7086270313451	25.7112772855539	25.7255982598915	25.7548926908193	25.8032519570040	25.8757908522466	25.9789714041452	26.1210578313829	26.3127691284544	26.5682373142862	26.9064526154258	27.3535111138434	27.9462379433694	28.7382800730748	29.8108847604867	31.2931868444701	33.4034805877761	36.5420777885591	41.5307782567411	50.3654813388633	69.4774823951482	137.454433620522	1436.73323747991	207.039567439610	97.0146963942635	56.9030239470375	36.7176155709928	24.9241746471234	17.4308680982016	12.4199984133465	8.96143912039753	6.53045877434873	4.80937510466554	3.59583994387590	2.75696450429921	2.20579659631698	1.89095159452867	1.79774625059558	1.96827814988314	2.57736821481699	4.27638775208220	11.2632348837092	42.4377543284410	9.99538290048535	6.77944598516005	5.72694904401283	5.30448548531993	5.14779393241938	5.12502276402959	5.17636609148263	5.27048623964025	5.38922908648720	5.52128440341820	5.65921854949713	5.79795103101160	5.93391242465023	6.06454975694628	6.18802071873476	6.30299615729760	6.40852757380479	6.50395523869941	6.58884259538730	6.66292821749981	6.72608981873983	6.77831674798037	6.81968859430930	6.85035828205845	6.87053852661190	6.88049084802116	6.88051656091521	6.87094931246867	6.85214884818649	6.82449576262473	6.78838704844011	6.74423229868320	6.69245044829651	6.63346696428027	6.56771141198628	6.49561533893011	6.41761042840424	6.33412688377119	6.24559201116169	6.15242897380028	6.05505569563041	5.95388389554064	5.84931823647531	5.74175557618029	5.63158430839042	5.51918378498993	5.40492381113425	5.28916420656114	5.17225442738111	5.05453324355362	4.93632846805301	4.81795673442623	4.69972332006218	4.58192201304193	4.46483502093261	4.34873292033410	4.23387464639459	4.12050752188324	4.00886732575171	3.89917840143207	3.79165380541071	3.68649549688505	3.58389456955337	3.48403152680466	3.38707660176426	3.29319012380680	3.20252293326591	3.11521684614529	3.03140517065904	2.95121327739427	2.87475922478571	2.80215444141267	2.73350446636223	2.66890974854195	2.60846650536266	2.55226764064143	2.50040372089487	2.45296400840497	2.41003754855015	2.40276145224860	2.44635108188700	2.49401782127648	2.54567528101941	2.60123013543417	2.66058119457866	2.72361855078943	2.79022281272208	2.86026443887152	2.93360318137161	3.01008764953960	3.08955500115792	3.17183076789214	3.25672881954749	3.34405147008756	3.43358972649955	3.52512367971424	3.61842303490198	3.71324777659597	3.80934896226991	3.90646963624984	4.00434585419861	4.10270780690761	4.20128103079132	4.29978769133288	4.39794792479554	4.49548122281513	4.59210784403562	4.68755023675267	4.78153445659137	4.87379156356073	4.96405898339012	5.05208181884640	5.13761409773533	5.22041994548198	5.30027467153376	5.37696576030522	5.45029375895189	5.52007305588803	5.58613254561504	5.64831617707075	5.70648338431413	5.76050939989602	5.81028545270973	5.85571885344413	5.89673297195887	5.93326711195398	5.96527628920532	5.99273092037789	6.01561643001220	6.03393278370419	6.04769395577611	6.05692733987027	6.06167311090117	6.06198354668721	6.05792231736355	6.04956375036764	6.03699207840237	6.02030067733406	5.99959130048629	5.97497331526037	5.94656294746047	5.91448253813781	5.87885981720352	5.83982719750338	5.79752109250559	5.75208126023342	5.70365017558084	5.65237243268679	5.59839417861416	5.54186257918526	5.48292531746775	5.42173012508329	5.35842434622586	5.29315453402650	5.22606607868482	5.15730286660391	5.08700696961077	5.01531836321965	4.94237467279510	4.86831094639626	4.79325945302889	4.71734950499619	4.64070730302115	4.56345580280936	4.48571460173051	4.40759984431709	4.32922414530850	4.25069652900610	4.17212238374833	4.09360343036407	4.01523770351444	3.93711954488928	3.85933960728128	3.78198486861963	3.70513865510358	3.62888067263549	3.55328704581059	3.47843036377800	3.40437973234338	3.40007672891801	3.44443810570082	3.49218787294088	3.54355424411085	3.59879080763900	3.65817947341525	3.72203385631818	3.79070317297735	3.86457674355100	3.94408920952428	4.02972660239738	4.12203342791641	4.22162096788017	4.32917704874892	4.44547758623682	4.57140029174176	4.70794102520134	4.85623340704870	5.01757246940600	5.19344334742446	5.38555630530724	5.59588978576862	5.82674370629742	6.08080595845396	6.36123608273569	6.67177151857703	7.01686385991296	7.40185547934038	7.83321118622678	8.31882600587984	8.86843993540733	9.49420570426953	10.2114796828179	11.0399454076055	12.0052451907327	13.1414096917298	14.4945811146996	16.1289118282632	18.1362812608607	20.6530632921879	23.8907335249412	28.1957696048328	34.1777595974500	43.0175556196253	57.3409839449379	84.4027262320209	154.321533982576	744.733060260958	281.076778227133	121.700828835871	79.1691217790383	59.5290463815154	48.2673795351555	41.0010202416942	35.9517840178773	32.2616389738366	29.4657503633183	27.2904455769869	25.5641447358837	24.1738480277670	23.0421305012640	22.1141895512029	21.3501681047065	20.7204058884622	20.2023936961492	19.7787583224469	19.4358929866273	19.1630042417433	18.9514347379747	18.7941729809645	18.6854924995909	18.6206822505925	18.5958424356783	18.6077279358642	18.6536268935627	18.7312655702535	18.8387330779991	18.9744213056573	19.1369765787496	19.3252604645406	19.5383177665492	19.7753502165628	20.0356947158819	20.3188052346095	20.6242376718980	20.9516371278659	21.3007271513470	21.6713006153887	22.0632119407820	22.4763704415206	22.9107346084151	23.3663071807036	23.8431308823495	24.3412847212888	24.8608807672971	25.4020613382843	25.9649965363382	26.5498820842814	27.1569374212627	27.7864040223197	30.4600367092051	34.1349011171074	38.4918564461062	43.7262161678744	50.1159703195285	58.0706693545001	68.2196399420931	81.5821712727546	99.9248172250107	126.599682198351	168.837583103628	245.640792785372	428.175314822525	1410.31176379477	1226.81190360654	445.749389151626	279.195052045565	206.801538503235	166.389221196865	140.650575293374	122.859275089284	109.856145546766	99.9613685702544	92.1993739382792	85.9648492063146	80.8623583213264	76.6226802665153	73.0562128750937	70.0256010209220	67.4289311544614	65.1890211106759	63.2463730413716	61.5544066329790	60.0761558450836	58.7819304773559	57.6476290631660	56.6535008433824	55.7832232982331	55.0232052378793	54.3620536399102	53.7901610518481	53.2993829190649	52.8827827860329	52.5344292902440	52.2492330799405	52.0228147968899	51.8513974430444	51.7317180432627	51.6609546944722	51.6366659715042	51.6567403230117	51.7193535950603	51.8229332064514	51.9661277983397	52.1477814129513	52.3669114381439	52.6226896979990	52.9144261834380	53.2415550076686	53.6036222441358	54.0002753634373	54.4312540333233	54.8963820847222	55.3955604785212	55.9287611339431	56.4960215009338	57.0974397768462	57.7331706825810	58.4034217257674	59.1084498889818	59.8485586897594	60.6240955665439	61.4354495509854	62.2830491923036	63.1673607039643	64.0888863067857	65.0481627459013	66.0457599618542	67.0822798985447	68.1583554328746	69.2746494127540	70.4318537917341	71.6306888499025	72.8719024918838	74.1562696138391	75.4845915322746	76.8576954682785	78.2764340815074	79.7416850488723	81.2543506834161	82.8153575893637	84.4256563497535	86.0862212434291	87.7980499885260	89.5621635098615	91.3796057279216	93.2514433673656	95.1787657831796	97.1626848028045	100];
% kappa = [22.0454076850486	21.9365958386363	21.8309609598464	21.7283615635978	21.6286579083878	21.5317119663279	21.4373873961616	21.3455495193167	21.2560652990413	21.1688033226761	21.0836337871065	21.0004284874452	20.9190608089874	20.8394057224855	20.7613397827847	20.6847411308646	20.6094894993278	20.5354662213793	20.4625542433386	20.3906381407272	20.3196041379763	20.2493401317985	20.1797357182692	20.1106822236664	20.0420727391179	19.9738021591075	19.9057672238972	19.8378665659214	19.7700007602186	19.7020723789644	19.6339860501811	19.5656485207002	19.4969687234634	19.4278578492537	19.3582294229584	19.2879993844728	19.2170861743669	19.1454108244459	19.0728970533510	18.9994713673612	18.9250631665721	18.8496048566474	18.7730319663594	18.6952832711552	18.6163009230155	18.5360305868973	18.4544215840857	18.3714270428177	18.2870040565778	18.2011138505142	18.1137219564731	18.0247983972091	17.9343178803899	17.8422600030945	17.7486094675788	17.6533563091837	17.5564961373648	17.4580303909431	17.3579666088184	17.2563187175406	17.1531073373187	17.0483601082499	16.9421120387904	16.8344058787607	16.7252925194908	16.6148314240717	16.5030910910960	16.3901495557539	16.2760949327114	16.1610260058469	16.0450528706834	15.9282976362407	15.8108951940707	15.6929940634571	15.5747573231998	15.4563636420970	15.3380084222439	15.2199050716499	15.1022864255088	14.9854063388421	14.8695414772928	14.7549933377198	14.6420905361257	14.5311914075693	14.4226869713675	14.3170043254481	14.2146105466496	14.1160171896807	14.0217854971321	13.9325324573729	13.8489378776819	13.7717526782707	13.7018086612005	13.6400300695967	13.5874473310091	13.5452134796943	13.5146238833333	13.4971400703654	13.4944186786810	13.5083468444431	13.4460907389537	13.0221260751535	12.6147249661972	12.2233560713230	11.8474993073696	11.4866458353724	11.1402980383273	10.8079694899467	10.4891849142668	10.1834801360017	9.89040202158330	9.60950841086395	9.34036803950886	9.08256045214958	8.83567590641844	8.59931526803563	8.37308989716933	8.15662152634169	7.94954213020395	7.75149378755417	7.56212853602001	7.38110821987632	7.20810433151125	7.04279784709717	6.88487905705961	6.73404739197193	6.59001124453220	6.45248778830280	6.32120279391156	6.19589044342594	6.07629314361753	5.96216133883468	5.85325332419414	5.74933505978934	5.65017998659415	5.55556884471460	5.46528949460982	5.37913674186626	5.29691216606685	5.21842395424963	5.14348673939975	5.07192144436390	5.00355513151913	4.93822085846889	4.87575753997820	4.81600981629893	4.75882792797483	4.70406759715548	4.65158991538979	4.60126123781249	4.55295308358309	4.50654204238579	4.46190968675148	4.41894248992001	4.37753174892161	4.33757351252287	4.29896851365290	4.26162210590156	4.22544420366232	4.19034922547851	4.15625604014287	4.12308791509647	4.09077246667408	4.05924161174875	4.02843152033862	3.99828256875307	3.96873929287387	3.93975034118854	3.91126842721862	3.88325028101322	3.85565659940940	3.82845199479360	3.80160494213335	3.77508772408502	3.74887637402127	3.72295061686018	3.69729380761792	3.67189286764579	3.64673821855279	3.62182371385382	3.59714656842268	3.57270728586733	3.54850958398181	3.52456031846551	3.50086940513494	3.47744974088633	3.45431712369869	3.43149017199643	3.40899024371826	3.38684135546410	3.36507010211502	3.34370557734161	3.32277929543463	3.30232511490673	3.28237916432749	3.26297977086345	3.24416739200212	3.22598455094317	3.20847577614061	3.19168754547819	3.20998915225446	3.28330052037588	3.36130587188711	3.44399639857752	3.53136061556160	3.62338356683316	3.72004601979295	3.82132364748227	3.92718619698079	4.03759664218733	4.15251031900301	4.27187404078124	4.39562519180310	4.52369079648482	4.65598656203388	4.79241589234719	4.93286887109928	5.07722121221014	5.22533317622321	5.37704845157900	5.53219300035487	5.69057386877613	5.85197796370780	6.01617079743658	6.18289520437196	6.35187003486505	6.52278883319288	6.69531850891640	6.86909801332448	7.04373703555527	7.21881473626826	7.39387854045122	7.56844301509985	7.74198886211289	7.91396206178970	8.08377320777237	8.25079708008206	8.41437250897090	8.57380258851198	8.72835530500283	8.87726465111778	9.01973230200942	9.15492993384548	9.28200226811652	9.40007092592328	9.50823917473266	9.60559764509918	9.69123108585792	9.76422621256598	9.82368068479758	9.86871322266840	9.89847484123030	9.91216114294832	9.90902556350907	9.88839341533467	9.84967651755989	9.79238814267879	9.71615795106219	9.62074652824381	9.50605909003041	9.37215788131248	9.21927277035666	9.04780953564492	8.85835536076618	8.65168109728276	8.42873992728751	8.19066215615004	7.93874598923851	7.67444428954195	7.39934746936500	7.11516283015867	6.82369082048739	6.52679882313078	6.22639319891822	5.92439039903100	5.62268800331797	5.32313654657786	5.02751295777929	4.73749636190211	4.45464688630363	4.18038798116519	3.91599261588529	3.66257355989752	3.42107780647828	3.19228505963922	2.97681008335365	2.77510861320635	2.58748645492651	2.41411134189256	2.25502709250705	2.11016959488246	1.97938414627509	1.86244368365903	1.75906745566540	1.66893970178154	1.59172792074488	1.52710032694045	1.47474211376136	1.43437017034381	1.40574593756941	1.40220622554357	1.41385730780322	1.43026356726339	1.45131805913290	1.47691406814061	1.50694556620691	1.54130770445862	1.57989734633874	1.62261364231303	1.66935864064469	1.72003792362016	1.77456125494000	1.83284322194790	1.89480385593921	1.96036921476245	2.02947191399476	2.10205159576957	2.17805532751143	2.25743792607362	2.34016220583334	2.42619915200119	2.51552802264116	2.60813638462923	2.70402009000655	2.80318319994367	2.90563786388392	3.01140416144779	3.12050991442778	3.23299047575428	3.34888850173075	3.46825371317418	3.59114265039810	3.71761842627509	3.84775048093913	3.98161434105354	4.11929138598855	4.26086862273062	4.40643847088614	4.55609855874342	4.70995153101787	4.86810486862003	5.03067072055176	5.19776574784546	5.36951097931005	5.54603167873109	5.72745722308414	5.91392099125748	6.10556026273825	6.30251612568957	6.50493339383520	6.71296053156685	6.92674958669742	7.14645613029755	7.37223920307261	7.60426126776015	7.84268816705314	8.08768908658164	8.33943652251324	8.59810625336031	8.86387731561038	9.13693198282246	9.41745574785896	9.70563730794759	10.0016685522916	10.3057445519700	10.6180635518898	10.9388269645738	11.2682393655841	11.6065084904012	11.9538452325937	12.3104636431275	12.6765809306800	13.0524174628341	13.4381967680424	13.8341455382584	14.2404936321455	14.6574740787809	15.0853230817807	15.5242800237801	15.9745874712093	16.4364911793118	16.9102400973580	17.3960863740116	17.8942853628107	18.4050956277311	18.9287789488021	19.4656003277478	20.0158279936327	20.5797334084896	21.1575912729140	21.7496795316087	22.3562793788660	22.9776752639760	23.6141548965509	24.2660092517591	24.9335325754598	25.6170223892339	26.3167794953069	27.0331079813598	27.7663152252252	30.4600367092051	34.1349011171074	38.4918564461062	43.7262161678744	50.1159703195285	58.0706693545001	68.2196399420931	81.5821712727546	99.9248172250107	126.599682198351	168.837583103628	245.640792785372	428.175314822525	1410.31176379477	1226.81190360654	445.749389151626	279.195052045565	206.801538503235	166.389221196865	140.650575293374	122.859275089284	109.856145546766	99.9613685702544	92.1993739382792	85.9648492063146	80.8623583213264	76.6226802665153	73.0562128750937	70.0256010209220	67.4289311544614	65.1890211106759	63.2463730413716	61.5544066329790	60.0761558450836	58.7819304773559	57.6476290631660	56.6535008433824	55.7832232982331	55.0232052378793	54.3620536399102	53.7901610518481	53.2993829190649	52.8827827860329	52.5344292902440	52.2492330799405	52.0228147968899	51.8513974430444	51.7317180432627	51.6609546944722	51.6366659715042	51.6567403230117	51.7193535950603	51.8229332064514	51.9661277983397	52.1477814129513	52.3669114381439	52.6226896979990	52.9144261834380	53.2415550076686	53.6036222441358	54.0002753634373	54.4312540333233	54.8963820847222	55.3955604785212	55.9287611339431	56.4960215009338	57.0974397768462	57.7331706825810	58.4034217257674	59.1084498889818	59.8485586897594	60.6240955665439	61.4354495509854	62.2830491923036	63.1673607039643	64.0888863067857	65.0481627459013	66.0457599618542	67.0822798985447	68.1583554328746	69.2746494127540	70.4318537917341	71.6306888499025	72.8719024918838	74.1562696138391	75.4845915322746	76.8576954682785	78.2764340815074	79.7416850488723	81.2543506834161	82.8153575893637	84.4256563497535	86.0862212434291	87.7980499885260	89.5621635098615	91.3796057279216	93.2514433673656	95.1787657831796	97.1626848028045	1000];
% kappa = [79.5863538049489	65.7278462494993	53.9170956049694	43.9145080404237	35.5001952405637	28.4731352789767	22.6503316077612	17.8659705439141	13.9705778246039	10.8301751038506	8.32543774169208	6.35085601496942	4.81390317277194	3.63421597148553	2.74279724156717	2.08125731548995	1.60112558575010	1.26329471668587	1.03773520506925	0.903821110029743	0.852231494990069	0.891641451743055	1.07374272604259	1.62196765712731	4.59374385652713	5.27781877837663	5.68498310476870	6.12780331095252	6.61144833552471	7.14250692262126	7.72947617127087	8.38346024317270	9.11919325880196	9.95657717609656	10.9230656725474	12.0574929333009	13.4164854413427	15.0857494039052	17.2011899208754	19.9915539060959	23.8734316278696	29.6929269773608	39.4697139663468	59.5078910197391	124.568136754709	952.724395066806	95.4108073455769	49.2206705753649	32.6776687112876	24.1694292970751	19.5312058977951	20.8716309469482	22.5799538338219	24.7768377291536	27.6466675317117	31.4853409570174	36.7987446616895	44.5280856087489	56.6448330071813	78.0886784235393	125.718823294506	318.893352142661	618.862886586464	159.454401438011	92.5475527052587	65.8067304070054	51.4773529329743	42.5872757230776	36.5621350919587	32.2304203600551	28.9826455531332	26.4705520866050	24.4808229669371	22.8756046670689	21.5619455011406	20.7288042170686	21.3679025927618	21.9156909116607	22.3647643283962	22.7087256970478	22.9423145784455	23.0615238178710	23.0637006081240	22.9476290289349	22.7135912282329	22.3634046900413	21.9004334172779	21.3295713408719	20.6571968404583	19.8910979108416	19.0403682139340	18.1152749944115	17.1271005835943	16.0879599430657	15.0105973816177	13.9081661931983	12.7939954913722	11.6813489471133	10.5831804705186	9.51189212572444	8.22284103311482	6.35704826673294	5.01000276127143	4.00738750345542	3.24458956064908	2.65500720658264	2.19421667536188	1.83145530339053	1.54476213086884	1.31806936970901	1.13939170769542	0.999661741023689	0.891960194993999	0.810994791785498	0.752739306275999	0.714177289525130	0.693114596560456	0.688037223401675	0.697999198625968	0.722530992865168	0.761562867227245	0.815360154894440	0.884468934755943	0.969671208958558	1.07194885130422	1.21092910909186	1.43284229947868	1.68963447395926	1.98459433462067	2.32150850158069	2.70485602535186	3.14007880871829	3.63396033553297	4.19516441493344	4.83501788473808	5.56867715650747	6.41691953679923	7.40899108578293	8.58732190669232	10.0157192591031	11.7944595414037	14.0901753801523	17.2007823439302	21.7150404734936	28.9793638312478	42.8884271071253	81.3050152027805	779.573259123744	100.930203094039	46.6137490594758	34.6445059756113	39.2235050158704	45.5087857727709	54.6798560228874	69.3273554980842	96.4637732806042	163.979745473999	624.536202418548	317.980489146563	122.522265753155	74.2183128522125	52.3358847024181	39.8536760909173	31.7870500334161	26.1467305611508	21.9827082178516	18.7839572201696	16.2511346998801	14.1972587326027	12.4995108262929	11.0738871559079	9.86098680496614	8.81761793278658	7.91162067870371	7.11855494815601	6.88360301808948	7.69512599348155	8.84803451455102	10.5980288208239	13.5422918736283	19.4751637044103	37.4015925613520	7752.47367400067	34.0195714157074	16.3372264399941	10.4463897582397	7.51309870624611	5.76407633591090	4.60733492662306	3.78895205410395	3.18190111686644	2.71558959386965	2.34765491863515	2.05113800307053	1.80806755399736	1.60600443897567	1.43606656356870	1.29174356707373	1.16815589092356	1.06157529267287	0.994488659769613	0.972071790793006	0.953516691375599	0.938755621255173	0.927752823352372	0.920505880188611	0.917047654191421	0.917448905180909	0.921821706983789	0.930323821959735	0.943164240163932	0.960610152987375	0.982995715150687	1.01073306207161	1.04432620371635	1.08438862842372	1.13166574660489	1.18706372349842	1.25168685180031	1.32688649177811	1.41432590590959	1.51606727651463	1.63469021567974	1.77345583354878	1.93653810798812	2.13110806810328	2.35143908800544	2.59483278248266	2.86604468114355	3.17151858955873	3.52015854926459	3.92456287734025	4.40307970303012	4.98340287319424	5.70924875924666	6.65370880768874	7.94862530685264	9.85803635509984	12.9985442905833	19.2210353878489	37.7867325799678	48646.3630849717	36.0964358114269	17.5836964141353	11.3869116298794	8.27405854298374	6.39772969418846	5.14204811184533	4.24308247939736	3.56888875598998	3.44663092041007	3.85622469359119	4.45762211328032	5.37774462189208	6.89467366298132	9.75536000698650	16.8768295726494	62.0591846912884	38.2232810355837	14.9390597442942	9.46062157989950	7.03950688328823	5.69094089482818	4.84251512203640	4.26778675705501	3.85933355609324	3.55973898689170	3.33561466677584	3.16629828334384	3.03838787443326	2.94287568876853	2.87354577713549	2.82603035686544	2.79722981025673	2.78494324199911	2.83870157333849	2.92911839163428	3.00809513919298	3.07509574206203	3.12982703224284	3.17222433865652	3.20243235062530	3.22078269412412	3.22776954955207	3.22402447477532	3.21029140383829	3.18740259040774	3.15625607120997	3.11779504887507	3.07298944183721	3.02281972379932	2.96826307670065	2.91028180726367	2.84981392507954	2.78776574639904	2.72500636873190	2.66236385364615	2.60062295572727	2.54052424183704	2.48276445437410	2.53756501224477	2.68087833906146	2.84457110871167	3.03124106797429	3.24410079406432	3.48715307481377	3.76543247896777	4.08534362111883	4.45514390677665	4.88564767500993	5.39127941116819	5.99169551387635	6.71436751266581	7.59886405714306	8.70429487809670	10.1230225328420	12.0078034797790	14.6306962157476	18.5276271327405	24.9204075518899	37.3296003719067	71.7794228936672	651.800971367806	95.9179459698924	45.5204475424284	40.3825308649332	43.2907285120606	46.3999026939634	49.7331700297562	53.3184685498158	57.1898352625698	61.3891021194339	65.9681798163898	70.9921821349814	76.5437758433627	82.7293568910519	89.6880155584149	97.6048804327469	106.731558795295	117.418506161482	130.168324542972	145.727685817166	165.255027199980	190.648498235169	225.246674522369	275.511918716991	355.794469849955	505.641922773984	888.342962859829	4012.60511077639	33.4162623379275	13.4077957735945	8.17231769702464	5.73193515741038	4.30432883492315	3.35853342696535	2.68100595609610	2.16946754196549	1.76892240307719	1.44724673582026	1.18445930144469	0.967487745846567	0.787401683738413	0.637852858419216	0.514149018981225	0.412681536624418	0.330562817543226	0.265397564156971	0.215148696090385	0.178080192201226	0.152771456386263	0.138201370570396	0.133892081041366	0.140083074234417	0.157887960338281	0.199486228188600	0.263199448138621	0.349581663397870	0.462832621221416	0.607836380460092	0.790416192983642	1.01777243739486	1.29918734972473	1.64716101452606	2.07929626900902	2.62156765363500	3.31432319332059	4.22412259332149	5.46933111238510	7.28262782304024	10.1922627218838	15.7073304720348	30.5549353776706	232.264887099839	44.5761705082622	20.7444741683869	13.5765756117335	10.0608878022472	7.93489998261034	6.48602637741386	6.39809842203567	6.59269225937882	6.79631166010926	7.01049232762885	7.23711627855651	7.47850475940833	7.73754273713566	8.01784827704039	8.32400695097483	8.66190245945014	9.03919300101706	9.46601438410059	9.95604675878715	10.5281851977955	11.2092546499800	12.0386198149582	13.0764363861821	14.4194086080552	16.2334443406075	18.8289489768278	22.8617815225084	29.9988410945611	46.0873674579784	116.281608390823	160.272463456525	824.916566976656	115.644300090719	50.5204831155907	30.9673732448296	21.6486297183311	16.2526572814310	12.7692821402689	10.3591154681152	8.60970538938383	7.29486105906427	6.28028618338337	5.48130176566559	4.84193972105496	4.32378367877793	3.89963769091008	3.54975013621578	3.25946577943123	3.01771340924112	2.81600176451234	2.64773521661584	2.50773658405575	2.39190762914401	2.29698319312384	2.22035033450901	2.15991343416793	2.12967868768994	2.10481214466119	2.08296292980621	2.06408234329397	2.04813762735480	2.03511071193443	2.02499719040830	2.01780549254690	2.01355622760589	2.01228167495852	2.01402540334420	2.01884200276320	2.02679691543630	2.03796635419460	2.05243729826100	2.07030755768740	2.09168589880361	2.11669222393312	2.14545779939933	2.17812552649786	2.21485025067862	2.25579910467238	2.30115188175283	2.35110143570956	2.40585410450130	1.84341552423356	1.51060621420170	1.30788540463703	1.17817547377568	1.09451164356845	1.04266740947318	1.01458635568410	1.00551777568238	1.01263059373773	1.03428742962518	1.06964011780954	1.11839292948707	1.18065834363263	1.25686633457765	1.34770585090127	1.45408632412933	1.57711201168933	1.71806477834321	1.87839255388485	2.05970168710818	2.26375202250708	2.49245390970193	2.74786660342248	3.03219767531104	3.0];
for i = 2:steps-1
    Vmax = 0.03;
    if(kappa(i)) < threshold
        Vmax = Vmax * (kappa(i) / threshold)^(1/2);
    end
    V = [V; Vmax];
end
V = [V; Vmax];

s_start = 0;
s_end = l_sum;
v_sup = V;
a_sup = 0.01;
j_max = 0.01;
% 加加速阶段
% t_jj = a_sup / j_max;
% v_jj = 1/2 * j_max * t_jj^2;
% s_jj = 1/6 * j_max * t_jj^3;
v0 = 0;
vt = 0;
V_reachability = velocity_planning(s_start, LL, v_sup, a_sup, v0, vt);
ssss = [s_start; LL];
a_reachability = 0;




% s_sum = [];
% v_sum = [];
% t_sum = [];

% s_sum_last = 0;
% v_sum_last = 0;
% t_sum_last = 0;
% a_m = [];
% V_real = 0;

% SS_L = 0;
% for i = 1:499
%     aa = (V_reachability(i+1)^2 - V_reachability(i)^2) /2 / L(i);
%     if aa > a_sup
%         a_m = [a_m, a_sup];
%     elseif aa < -a_sup
%         a_m = [a_m, -a_sup];
%     else
%         a_m = [a_m, aa];
%     end
%     delta_t = (V_reachability(i+1) - V_reachability(i)) / a_m(i);
%     vvv = min(delta_t* a_m(i)+V_reachability(i), V_reachability(i+1));
%     t1 = delta_t* a_m(i) + V_reachability(i);
%     t2 = V_reachability(i+1);
% 
%     V_real = [V_real, vvv];
% end
SSSS = 0;
sum_T = 0;
sum_S = 0;
TTTT = 0;
% sum_V = 0;
VVVV = 0;
AAAA = 0;
allow_omega = 1;
scale = 1;
for f = 1:10
SSSS = 0;
sum_T = 0;
sum_S = 0;
TTTT = 0;
VVVV = 0;
AAAA = 0;
V_reachability = V_reachability * scale;
Vmax = Vmax * scale;
for i=2:steps-1    
    if(V_reachability(i+1)-V_reachability(i)<0) && (V_reachability(i)-V_reachability(i-1)>=0)
        List = [0, 2*ssss(i), Vmax, a_sup, j_max];
        [Tp, Jp, Ap, Vp, Sp] = P2PMultiAxisDoubleSTrajectory(0.001, List);
        break;
    end
end

for i=steps-1:-1:2    
    if(V_reachability(i)-V_reachability(i-1)>0) && (V_reachability(i+1)-V_reachability(i)<=0)
        List = [0, 2*(ssss(end)-ssss(i)), Vmax, a_sup, j_max];
        [Tpf, Jpf, Apf, Vpf, Spf] = P2PMultiAxisDoubleSTrajectory(0.001, List);
        break;
    end
end
% s_type_ssss = 0;
% for jj = 1:steps-1
%     if ssss(jj) < s_jj && ssss(jj+1) > s_jj
%         s_type_ssss = [s_type_ssss, s_jj]
%     elseif ssss(jj) > s_jj 
%         s_type_ssss = [s_type_ssss, ssss(jj)];
%     end
% end
% s_type_ssss = [s_type_ssss, ssss(steps)];

for i = 1: steps-1
v0 = V_reachability(i);
vt = V_reachability(i+1);
v_max = max(V_reachability(i), V_reachability(i+1));
% a_max = 0.3;
s = ssss(i+1) - ssss(i);
% [s_i, v_i, t_i] = trapezoidal_velocity_profile(v0, vt, v_max, a_max, s);
% 
% s_sum = [s_sum, s_sum_last + s_i];
% v_sum = [v_sum, v_i];
% t_sum = [t_sum, t_sum_last + t_i];
% 
% s_sum_last = s_sum(end);
% v_sum_last = v_sum(end);
% t_sum_last = t_sum(end);
% SS_L = [SS_L; s_sum_last];


TTT = 2 * s  / (v0+ vt);
sum_T = sum_T + TTT;
a = (vt - v0) / TTT;
s = v0 * TTT + 0.5 * a * TTT^2;
% if (a > a_sup)
%     a = a_sup;
%     vt = v0 + TTT * a;
% %     t1 = TTT
%     s = v0 * TTT + 0.5 * a * TTT^2;
%     sum_S = sum_S + s;
% SSSS = [SSSS, sum_S];
% TTTT = [TTTT, sum_T];
% VVVV = [VVVV, vt];
% AAAA = [AAAA, a];
% elseif (a  < -a_sup)
%     a = -a_sup;
%     s = s + SSSS(i-1) - SSSS(i-2);
%     v0 = VVVV(i-1)
%     TTT = 2 * ss  / (v0+ vt);
%     sum_T = sum_T + TTT;
%     a = (vt - v0) / TTT;
%     s = v0 * TTT - 0.5 * a * TTT^2;
% 
%     sum_S = sum_S + s;
%     SSSS = [SSSS, sum_S];
%     TTTT = [TTTT, sum_T];
%     VVVV = [VVVV, vt];
%     AAAA = [AAAA, a];
% else


%  % 先梯形规划，如果加速度超标，那就加入超标的前一段和后一段重新规划
%     sum_S = sum_S + s;
% SSSS = [SSSS, sum_S];
% TTTT = [TTTT, sum_T];
% VVVV = [VVVV, vt];
% AAAA = [AAAA, a];
% end

    sum_S = sum_S + s;
SSSS = [SSSS, sum_S];
TTTT = [TTTT, sum_T];
VVVV = [VVVV, vt];
AAAA = [AAAA, a];

end

%% S型速度加速
[~,w] = size(Ap);
for jj=2:w
    if abs(Jp(jj)) < 1e-6 && Jp(jj-1) < 0
        TTp = Tp(1:jj);
        SSp = Sp(1:jj);
        VVp = Vp(1:jj);
        AAp = Ap(1:jj);
        break;
    end
end
b1 = [];
deltaT = 0;
NTT = [];
for jj =2:steps
    deltaT = [deltaT, TTTT(jj)-TTTT(jj-1)];
end
for jj =2:steps
    if SSSS(jj)>SSp(end) && SSSS(jj-1)<SSp(end)
        deltaS = SSSS(jj) - SSp(end);
        Tjj = 2 * deltaS / (VVVV(jj) + VVp(end));
        for ii=jj:steps
            deltaT(jj) = Tjj;
            NTT = [NTT, TTp(end) + sum(deltaT(jj:ii))];
        end
        SSSS = [SSp, SSSS(jj:end)];
        VVVV = [VVp, VVVV(jj:end)];
        AAAA = [AAp, AAAA(jj:end)];
        TTTT = [TTp, NTT];
        begins = jj;
        break;
    end
end

%% S型速度减速 
[~,llen] = size(AAAA);
[~,w] = size(Apf)
for jj=w-1:-1:2
    if abs(Jpf(jj)) < 1e-6 && Jpf(jj+1) < 0
        TTpf = Tpf(jj:end);
        SSpf = Spf(jj:end);
        VVpf = Vpf(jj:end);
        AApf = Apf(jj:end);
        break;
    end
end
deltas = SSpf(end) - SSpf(1);
[~,lenpf] = size(AApf);
NTT = [];
NSS = [];
for jj =llen:-1:2
    if (SSSS(end)-SSSS(jj))<deltas && (SSSS(end)-SSSS(jj-1))>deltas
        deltaS = SSSS(end) - SSSS(jj-1) - deltas;
        deltaT = deltaS / VVpf(1);
        importS = SSSS(jj-1) + deltaS;
        importT = TTTT(jj-1) + deltaT;
        Tjj = 2 * deltaS / (VVVV(jj) + VVp(end));
        for ii=1:lenpf
            NSS = [NSS, importS + SSpf(ii) - SSpf(1)];
            NTT = [NTT, importT + TTpf(ii) - TTpf(1)];
        end
        SSSS = [SSSS(1:jj-1), importS, NSS(2:end)];
        VVVV = [VVVV(1:jj-1), VVpf(1), VVpf(2:end)];
        AAAA = [AAAA(1:jj-1), AApf(1), AApf(2:end)];
        TTTT = [TTTT(1:jj-1), importT, NTT(2:end)];
        ends = llen - jj;
        break;
    end
end

%% 不同规划方式融合
[mm,nn] = size(TTTT);
Cs = [];
uus = [];
dt = 0.05;
for i = 0:dt:TTTT(end)
    for j = 1:length(TTTT)-1
        if (i >= TTTT(j)) && (i < TTTT(j+1))

            delta = i - TTTT(j);
            percent = delta/(TTTT(j+1) - TTTT(j));
            lucheng = percent * (SSSS(j+1) - SSSS(j)) + SSSS(j);
            
            uu = lucheng / SSSS(end);
            if(uu > 1)
            uu = 1;
            end
        end
    end

    uus = [uus, uu];
end

for i = 1:1:length(uus)
    if(uus(i) == 1)
        uus(i) = uus(i) - 1e-12;
    end
    [N0, N1, N2, N3] = bspline_basis(p, knots, uus(i));
    Css = N0 * P(:,1:3);
    Cs = [Cs; Css];
end

  q= [];
for i = 1:n+1
qq = eul2quat(P(i,4:6), "ZYX");
q = [q; qq];
end
num_samples = length(uus);
t = linspace(0, 1, num_samples);
k = 1;
s = sCurve(t, k);
VAL = [q(1,:)];
for i =2:num_samples-1
    val = quat_squad(q',s(i));
    VAL = [VAL; val];
end
VAL = [VAL; q(n+1,:)];
dT = ones(num_samples,1);
dT = dT * dt;
% dt(1) = 0.05;
% dt(1000) = 0.05;
T = 0;
[omega, alpha] = computeAngularVelocityAndAcceleration(VAL, dT);
for i = 1:num_samples-1
    T = [T,sum(dT(1:i))];
end


nomal_omega = [];
for i = 1:length(omega)
    nomega = sqrt(omega(i,1)^2 + omega(i,2)^2 + omega(i,3)^2);
    nomal_omega = [nomal_omega, nomega];
end
id=find(nomal_omega==max(nomal_omega));
omega_max=nomal_omega(id);
scale = allow_omega / omega_max;
% plot(T, nomal_omega(:));
% hold on;
% plot(T, omega(:,2));
% hold on;
% plot(T, omega(:,3));
EULER = [];
for i = 1:length(VAL)
    euler = quat2eul(VAL(i,:), "ZYX");
    EULER = [EULER; euler];
end

if scale >= 1
    break;
end
end
za = [0,0,1]';
data4 = [];
for i = 1:num_samples
    R = quat2rotm(VAL(i,:));
    data4 = [data4; (R*za)'];
end
%      plot3(data1(:, 1), data1(:, 2), data1(:, 3), 'b', 'LineWidth', 1.5);
%     hold on;
%      plot3(data2(:, 1), data2(:, 2),data2(:, 3), 'ro--', 'LineWidth', 0.7, 'MarkerSize', 4);
%      legend('B 样条曲线', '控制点');
%     grid on;
%     hold on;

figure;
plot3(C(:, 1), C(:, 2), C(:, 3), 'b', 'LineWidth', 1.5);
hold on;
plot3(P(:, 1), P(:, 2),P(:, 3), 'ro--', 'LineWidth', 1, 'MarkerSize', 8);
grid on;
hold on;
plot3(Cs(:, 1), Cs(:, 2),Cs(:, 3), 'go', 'LineWidth', 1, 'MarkerSize', 8);
quiver3(Cs(:, 1), Cs(:, 2),Cs(:, 3), data4(:,1),data4(:,2),data4(:,3));

figure;
subplot(3,1,1);
plot(T, VAL(:,1));
hold on;    
plot(T, VAL(:,2));
plot(T, VAL(:,3));
plot(T, VAL(:,4));
legend('加速度变化');
subplot(3,1,2);
plot(T, omega(:,1));
hold on;    
plot(T, omega(:,2));
plot(T, omega(:,3));
legend('加速度变化');
subplot(3,1,3);
plot(T, alpha(:,1));
hold on;    
plot(T, alpha(:,2));
plot(T, alpha(:,3));
legend('角加速度变化');

figure;
subplot(3,1,1);
plot(TTTT,  SSSS, 'LineWidth', 2);
legend('位置变化');
subplot(3,1,2);
plot(TTTT,  VVVV, 'LineWidth', 2);
legend('速度变化');
subplot(3,1,3);
plot(TTTT,  AAAA, 'LineWidth', 2);
legend('迹连续变化');


figure;
subplot(3,1,1);
plot(T,  EULER(:,1), 'LineWidth', 2);
legend('x');
subplot(3,1,2);
plot(T,  EULER(:,2), 'LineWidth', 2);
legend('y');
subplot(3,1,3);
plot(T,  EULER(:,3), 'LineWidth', 2);
legend('z');
