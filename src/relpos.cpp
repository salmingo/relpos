/*
 Name        : relpos.cpp
 Author      : Xiaomeng Lu
 Version     :
 Copyright   : SVOM Group, NAOC
 Description : 计算GWAC系统中, JFoV相对FFoV的位置
 1) 输入参数: 文件名1 文件名2 经度方向基准 倾斜方向基准
    文件内容分为三列, 分别对应: 赤经 赤纬 FITS文件名
    赤经/赤纬: 量纲: 角度
    文 件 名:  G<cam_id>_<imgtypabbr>_<utc>.fit
              cam_id : 相机标志, 三字节字符串
              imgtypabbr: 图像类型缩写
              utc    : UTC时间, 格式: YYMMDDThhmmssfs, fs量纲为10毫秒
    方向基准量纲: 角度. 缺省时为0
 2) cam_id整数为5倍数对应FFoV
 3) 查找与JFoV时间点最接近FFoV的位置
 4) 计算各时间点JFoV相对FFoV的位置
 5) 输出: JFoV原始位置, JFoV文件名, FFoV参考位置, FFoV文件名, 相对经度, 相对倾斜, 相对基准经度, 相对基准倾斜到控制台
    输出: JFoV原始位置, JFoV文件名, FFoV参考位置, FFoV文件名, 相对经度, 相对倾斜, 相对基准经度, 相对基准倾斜到文件
    输出文件名格式: G<cam_id>_<hhmm>-<hhmm>.txt
                  cam_id为JFoV相机标志
                  第一个hhmm为JFoV起始时间
                  第二个hhmm为JFoV结束时间

 约束条件:
 1) 望远镜处于跟踪模式
 2) 文件数据从前到后按照时间顺序
 */

#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

using std::string;
using std::vector;

//////////////////////////////////////////////////////////////////////////////
/// 宏定义
#define API		3.141592653589793
#define PI360	6.283185307179586
#define D2R		0.017453292519943		// 使用乘法, 角度转换为弧度的系数
#define R2D		57.295779513082323		// 使用乘法, 弧度转换为角度的系数
#define reduce(x, period)	((x) - floor((x) / (period)) * (period))

void Sphere2Cart(double r, double alpha, double beta, double& x, double& y, double& z)
{
	x = r * cos(beta) * cos(alpha);
	y = r * cos(beta) * sin(alpha);
	z = r * sin(beta);
}

void Cart2Sphere(double x, double y, double z, double& r, double& alpha, double& beta)
{
	r = sqrt(x * x + y * y + z * z);
	if ((alpha = atan2(y, x)) < 0) alpha += PI360;
	beta  = atan2(z, sqrt(x * x + y * y));
}

void RotateForward(double alpha0, double beta0, double& alpha, double& beta)
{
	double r = 1.0;
	double x1, y1, z1;	// 原坐标系投影位置
	double x2, y2, z2;	// 新坐标系投影位置

	// 在原坐标系的球坐标转换为直角坐标
	Sphere2Cart(r, alpha, beta, x1, y1, z1);
	/*! 对直角坐标做旋转变换. 定义矢量V=(alpha0, beta0)
	 * 主动视角, 旋转矢量V
	 * 先绕Z轴逆时针旋转: -alpha0, 将矢量V旋转至XZ平面
	 * 再绕Y轴逆时针旋转: -(PI90 - beta0), 将矢量V旋转至与Z轴重合
	 **/
	 x2 = sin(beta0) * cos(alpha0) * x1 + sin(beta0) * sin(alpha0) * y1 - cos(beta0) * z1;
	 y2 = -sin(alpha0) * x1 + cos(alpha0) * y1;
	 z2 = cos(beta0) * cos(alpha0) * x1 + cos(beta0) * sin(alpha0) * y1 + sin(beta0) * z1;
	// 将旋转变换后的直角坐标转换为球坐标, 即以(alpha0, beta0)为极轴的新球坐标系中的位置
	Cart2Sphere(x2, y2, z2, r, alpha, beta);
}
//////////////////////////////////////////////////////////////////////////////
/// 数据结构
struct PointRaw {// 原始单数据点
	double ra, dc;	//< 赤经, 赤纬, 量纲: 角度
	int ymd;			//< 年月日
	int hh, mm, ss;	//< 时分秒, 秒量纲: 0.01秒
	double secs;		//< 秒数
	string fname;	//< 文件名

public:
	PointRaw& operator=(const PointRaw& other) {
		if (this != &other) {
			ra = other.ra;
			dc = other.dc;
			ymd = other.ymd;
			hh = other.hh;
			mm = other.mm;
			ss = other.ss;
			fname = other.fname;
		}

		return *this;
	}
};
typedef vector<PointRaw> PtRV;	//< 原始数据点集合

struct PointFile {// 文件数据点
	string cid;		//< 相机标志
	PtRV pts;		//< 数据点集合

public:
	virtual ~PointFile() {
		pts.clear();
	}
};

struct PointCross {// 交叉数据点
	double ra, dc;	//< JFoV中心位置, 量纲: 角度
	string fname;	//< JFoV文件名
	double ra0, dc0;	//< FFoV中心位置, 量纲: 角度
	string fname0;	//< FFoV文件名
	double rot, tilt;	//< 旋转角和倾斜角, 量纲: 角度

public:
	/*!
	 * @brief 设置数据点, 即JFoV数据
	 * @param pt 原始数据
	 */
	void SetPoint(const PointRaw& pt) {
		ra = pt.ra;
		dc = pt.dc;
		fname = pt.fname;
	}

	/*!
	 * @brief 设置参考点, 即FFoV数据
	 * @param pt 原始数据
	 */
	void SetPointRef(const PointRaw& pt) {
		ra0 = pt.ra;
		dc0 = pt.dc;
		fname0 = pt.fname;

		rot = ra * D2R;
		tilt= dc * D2R;
		RotateForward(ra0 * D2R, dc0 * D2R, rot, tilt);
		rot *= R2D;
		tilt = 90 - tilt * R2D;
	}
};
//////////////////////////////////////////////////////////////////////////////
/// 全局变量
string pathSrc1, pathSrc2;	//< 输入文件路径名
double rot0, tilt0;	//< 旋转与倾斜基准, 量纲: 角度
PointFile pt_jfov, pt_ffov;	//< JFoV和FFoV的文件数据集合
bool bjfov, bffov;	//< JFoV和FFoV数据完备标志
string pathDst; //< 输出文件名
vector<PointCross> pt_cross;		//< 数据交叉结果

//////////////////////////////////////////////////////////////////////////////
/// 子函数
/*!
 * @brief 解析行信息
 * @param line  行信息
 * @param ra    赤经
 * @param dc    赤纬
 * @param fname 文件名
 */
void ResolveLine(const char* line, double& ra, double& dc, string& fname) {
	char* token;
	char seps[] = " \t\r\n";
	int n = strlen(line);
	char* buff = new char[n + 1];

	strcpy(buff, line);
	token = strtok(buff, seps);  ra = atof(token);
	token = strtok(NULL, seps);  dc = atof(token);
	token = strtok(NULL, seps);  fname = token;
	delete []buff;
}

/*!
 * @brief 解析文件名
 * @param fname  文件名
 * @param cid    相机标志
 * @param ymd    年月日
 * @param hms    时分秒
 */
void ResolveFilename(const char* fname, string& cid, int& ymd, int& hms) {
	char* token;
	char seps[] = "G_T.fit";
	int pos(0);
	int n = strlen(fname);
	char* buff = new char[n + 1];
	strcpy(buff, fname);

	token = strtok(buff, seps);
	while(token) {
		switch(++pos) {
		case 1: // cid
			cid = token;
			break;
		case 2: // obstyp or imgtyp
			if (strcasecmp(token, "mon")) ++pos;
			break;
		case 3: // imgtyp
			break;
		case 4: // ymd
			ymd = atoi(token);
			break;
		case 5: // hms
			hms = atoi(token);
			break;
		default:
			break;
		}

		token = strtok(NULL, seps);
	}

	delete []buff;
}

/*!
 * @brief 解析文件内容
 * @param filepath 原始文件路径
 * @return
 * 文件解析结果
 */
bool ResolveFile(const string& filepath) {
	printf("\n");

	FILE *fp = fopen(filepath.c_str(), "r");
	if (!fp) return false;
	printf("---------- Resolving file: %s ----------\n", filepath.c_str());

	char line[200];	// 行缓存区
	double ra, dc;	// 赤经/赤纬
	string fname;	// 文件名
	string cid;		// 相机标志
	int ymd, hms;	// 时间
	int hh, mm, ss;	// 时分秒
	int n(0);
	PointFile* ptr = NULL;	// 文件数据指针
	bool* valid;

	while(!feof(fp)) {
		if (!fgets(line, 200, fp)) continue;
		ResolveLine(line, ra, dc, fname);
		ResolveFilename(fname.c_str(), cid, ymd, hms);
		if (!ptr) {
			if (atoi(cid.c_str()) % 5 == 0) {
				ptr = &pt_ffov;
				valid = &bffov;
				printf("file<%s> is considered to be from FFoV\n", filepath.c_str());
			}
			else {
				ptr = &pt_jfov;
				valid = &bjfov;
				printf("file<%s> is considered to be from JFoV\n", filepath.c_str());
			}
			ptr->cid = cid;
		}
		ss = hms % 10000;
		hms /= 10000;
		mm = hms % 100;
		hh = hms / 100;
		++n;
//		printf("%8.4f %8.4f %s %s %06d %02d %02d %04d\n", ra, dc, fname.c_str(), cid.c_str(), ymd, hh, mm, ss);

		PointRaw pt;
		pt.ra = ra;
		pt.dc = dc;
		pt.fname = fname;
		pt.ymd = ymd;
		pt.hh = hh;
		pt.mm = mm;
		pt.ss = ss;
		pt.secs = (hh * 60 + mm) * 60 + ss * 0.01;
		ptr->pts.push_back(pt);
	}
	fclose(fp);

	printf("%d points are resolved from file\n", n);
	*valid = n > 0;
	if (ptr == &pt_jfov) {
		char buff[100];
		sprintf(buff, "G%s_%02d%02d-%02d%02d.txt", ptr->cid.c_str(),
				ptr->pts[0].hh, ptr->pts[0].mm,
				ptr->pts[n - 1].hh, ptr->pts[n - 1].mm);
		pathDst = buff;
	}

	return (n > 0);
}

/*!
 * @brief 检查JFoV或FFoV时间有效性
 * @return
 * 时间有效性
 * @note
 * 有效性判据: 数据日期相同
 */
bool TimeCheck(const PointFile* ptf) {
	int n = ptf->pts.size(), i;
	int ymd = ptf->pts[0].ymd;
	for (i = 1; i < n && ymd == ptf->pts[i].ymd; ++i);
	return (i == n);
}

/*!
 * @brief 检查JFoV和FFoV的时间有效性
 * @return
 * 时间有效性
 * @note
 * 有效性判据: JFoV和FFoV时间范围有交集
 * @note
 * 2017-10-28
 * 为了通用性, 应做更精细判读. 例如: 允许单点数据
 */
bool TimeCrossCheck() {
	return (pt_ffov.pts[0].ymd == pt_jfov.pts[0].ymd);
}

/*!
 * @brief 从FFoV原始数据中找到与秒数最接近的数据点
 * @param secs JFoV秒数
 * @param from 起始扫描位置
 * @param n    FFoV数据长度
 * @return
 * 匹配数据点位置
 *  -1: 未找到匹配数据
 * >=0: 匹配数据位置
 * @note
 * 匹配条件: 秒数相差不超过10
 */
int FindMatchedData(double secs, int from, int n) {
	PointRaw* pt;
	int i;
	double dt0(fabs(secs - pt_ffov.pts[from].secs));
	double dt1(1E30);

	for (i = from + 1; i < n; ++i) {
		dt1 = fabs(secs - pt_ffov.pts[i].secs);
		if (dt1 > dt0) break;
		dt0 = dt1;
		from = i;
	}

	return (fabs(secs - pt_ffov.pts[from].secs) > 10.0 ? -1 : from);
}

/*!
 * @brief 扫描原始数据并计算相对位置并输出结果
 */
void ScanData() {
	printf("\nscan and try to find matched data\n");

	int n1 = pt_jfov.pts.size();
	int n2 = pt_ffov.pts.size();
	int i, j(0), k;
	PointRaw* pt;

	for (i = 0; i < n1; ++i) {
		pt = &pt_jfov.pts[i];
		if ((k = FindMatchedData(pt->secs, j, n2)) >= 0) {
			j = k;

			PointCross ptc;
			ptc.SetPoint(*pt);
			ptc.SetPointRef(pt_ffov.pts[j]);
			pt_cross.push_back(ptc);
		}
	}

	printf("found %lu matched points\n", pt_cross.size());
}

/*!
 * @brief 输出处理结果到文件
 * @param fp 文件描述符
 * @note
 * fp==stdout或stderr时, 输出到控制台
 */
void OutputResult(FILE* fp) {
	printf("\n");
	int n = pt_cross.size(), i;
	PointCross* pt;
	double drot;
	if (n == 0) return;

	fprintf(fp, "%8s %8s %33s %8s %8s %33s %5s %4s %6s %5s\n",
			"R.A.  ", "DEC.  ", "FileName            ",
			"R.A.0 ", "DEC.0 ", "FileName.0          ",
			"Rot ", "Tilt", "rRot ", "rTilt");
	for (i = 0; i < n; ++i) {
		pt = &pt_cross[i];
		drot = rot0 - pt->rot;
		if (drot > 180.0) drot -= 360.0;
		else if (drot < -180.0) drot += 360.0;

		fprintf(fp, "%8.4f %8.4f %33s %8.4f %8.4f %33s %5.1f %4.1f %6.1f %5.1f\n",
				pt->ra, pt->dc, pt->fname.c_str(),
				pt->ra0, pt->dc0, pt->fname0.c_str(),
				pt->rot, pt->tilt, drot, tilt0 - pt->tilt);
	}

	if (fp == stdout || fp == stderr) {
		double rsum(0.0), rsq(0.0), tsum(0.0), tsq(0.0);
		double rmin(1E30), rmax(-1E30), tmin(1E30), tmax(-1E30);
		double rmean, rrms, tmean, trms;
		double rot(pt_cross[0].rot), tilt;

		for (i = 0; i < n; ++i) {
			pt = &pt_cross[i];
			tilt = pt->tilt;
			drot = pt->rot - rot;
			if (drot > 180.0) rot = pt->rot - 360.0;
			else if (drot < -180.0) rot = pt->rot + 360.0;
			else rot = pt->rot;

			if (rmin > rot) rmin = rot;
			if (rmax < rot) rmax = rot;
			if (tmin > tilt) tmin = tilt;
			if (tmax < tilt) tmax = tilt;

			rsum += rot;
			rsq += (rot * rot);
			tsum += tilt;
			tsq += (tilt * tilt);
		}

		rmean = rsum / n;
		tmean = tsum / n;
		rrms = sqrt((rsq - rsum * rmean) / n);
		trms = sqrt((tsq - tsum * tmean) / n);
		printf("****************************** Statistical results ******************************\n");
		printf("Rotation Minimum = %6.1f \t Rotation Maximum = %6.1f\n", reduce(rmin, 360.0), reduce(rmax, 360.0));
		printf("Rotation Mean    = %6.2f \t Rotation Stdev   = %6.2f\n", reduce(rmean, 360.0), rrms);
		printf("Tilt Minimum     = %6.1f \t Tilt Maximum     = %6.1f\n", tmin, tmax);
		printf("Tilt Mean        = %6.2f \t Tilt Stdev       = %6.2f\n", tmean, trms);
		printf("****************************** Statistical results ******************************\n");
	}
}
//////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
	if (argc < 3) {
		printf("\nUsage:\n\trelpos <path 1> <path 2> <rotation base> <inclination base>\n");
		return -1;
	}
	pathSrc1 = argv[1];
	pathSrc2 = argv[2];
	rot0  = argc >= 4 ? atof(argv[3]) : 0.0;
	tilt0 = argc >= 5 ? atof(argv[4]) : 0.0;
	bjfov = bffov = false;

	if (!ResolveFile(pathSrc1)) {
		printf("\nfail to resolve file<%s>\n", pathSrc1.c_str());
		return -2;
	}
	if (!ResolveFile(pathSrc2)) {
		printf("\nfail to resolve file<%s>\n", pathSrc2.c_str());
		return -2;
	}
	if (!bjfov) {
		printf("\nJFoV data is unavailable\n");
		return -3;
	}
	if (!bffov) {
		printf("\nFFoV data is unavailable\n");
		return -3;
	}
	if (!(TimeCheck(&pt_jfov) && TimeCheck(&pt_ffov) && TimeCrossCheck())) {
		printf("\ntime range do not match\n");
		return -4;
	}

	ScanData();

	if (!pt_cross.size()) {
		printf("\nno any data matches condition\n");
	}
	else {
		OutputResult(stdout); // 输出到控制台
		// 输出到文件
		FILE* fp = fopen(pathDst.c_str(), "w");
		if (fp) {
			printf("---------- results are saved as file<%s> ----------\n", pathDst.c_str());
			OutputResult(fp);
			fclose(fp);
		}
		else {
			printf("\nfailed to create result file<%s>\n", pathDst.c_str());
		}

		pt_cross.clear();
	}

	printf("\n");

	return 0;
}
