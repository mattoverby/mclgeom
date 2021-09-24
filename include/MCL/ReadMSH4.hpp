// From IPC: https://github.com/ipc-sim/IPC/blob/master/src/Utils/IglUtils.hpp
// License: MIT

#ifndef MCL_READMSH4_HPP
#define MCL_READMSH4_HPP 1

#include <fstream>
#include <Eigen/Core>

namespace mcl
{

// Reads a version 4 MSH file, not supported by igl.
static inline bool readMSH4(const std::string& filePath, Eigen::MatrixXd& TV, Eigen::MatrixXi& TT)
{
    FILE* in = fopen(filePath.c_str(), "r");
    if (!in) {
        return false;
    }

    TV.resize(0, 3);
    TT.resize(0, 4);
	char* tmp;
	int tmp_i;

    char buf[BUFSIZ];
    while ((!feof(in)) && fgets(buf, BUFSIZ, in)) {
        if (strncmp("$Nodes", buf, 6) == 0) {
            tmp = fgets(buf, BUFSIZ, in);
            int vAmt;
            sscanf(buf, "1 %d", &vAmt);
            TV.resize(vAmt, 3);
            tmp = fgets(buf, BUFSIZ, in);
            break;
        }
    }
    
    if (TV.rows()== 0)
        return false;

    int bypass;
    for (int vI = 0; vI < TV.rows(); vI++) {
        tmp_i = fscanf(in, "%d %le %le %le\n", &bypass, &TV(vI, 0), &TV(vI, 1), &TV(vI, 2));
    }

    while ((!feof(in)) && fgets(buf, BUFSIZ, in)) {
        if (strncmp("$Elements", buf, 9) == 0) {
            tmp = fgets(buf, BUFSIZ, in);
            int elemAmt;
            sscanf(buf, "1 %d", &elemAmt);
            TT.resize(elemAmt, 4);
            tmp = fgets(buf, BUFSIZ, in);
            break;
        }
    }

    if (TT.rows() == 0)
        return false;

    int minTT = 9999;
    for (int elemI = 0; elemI < TT.rows(); elemI++) {
        tmp_i = fscanf(in, "%d %d %d %d %d\n", &bypass,
            &TT(elemI, 0), &TT(elemI, 1), &TT(elemI, 2), &TT(elemI, 3));
        minTT = std::min(minTT, TT.row(elemI).minCoeff());
    }
    
    if (minTT != 0) {
        TT.array() -= minTT;
    }

	(void)(tmp);
	(void)(tmp_i);

    // while ((!feof(in)) && fgets(buf, BUFSIZ, in)) {
    //     if (strncmp("$Surface", buf, 7) == 0) {
    //         fgets(buf, BUFSIZ, in);
    //         int elemAmt;
    //         sscanf(buf, "%d", &elemAmt);
    //         F.resize(elemAmt, 3);
    //         break;
    //     }
    // }
    // for (int triI = 0; triI < F.rows(); triI++) {
    //     fscanf(in, "%d %d %d\n", &F(triI, 0), &F(triI, 1), &F(triI, 2));
    // }

    fclose(in);

    return TV.rows() > 0 && TT.rows() > 0;
}

} // namespace mcl

#endif
