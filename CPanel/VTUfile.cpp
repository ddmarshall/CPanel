//
//  VTUfile.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 11/25/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "VTUfile.h"

void VTUfile::write()
{
    std::ofstream f;
    f.open(name);
    if (f.is_open())
    {
        f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        f << "\t<UnstructuredGrid>\n";
        for (pieces_index_type p=0; p<pieces.size(); p++)
        {
            piece piece = pieces[p];
            f << "\t\t<Piece NumberOfPoints=\"" << piece.pnts.rows() << "\" NumberOfCells=\"" << piece.connectivity.rows() << "\">\n";
            if (piece.cellData.size() != 0)
            {
                f << "\t\t\t<CellData Scalars=\"scalars\">\n";
                for (pieces_index_type c=0; c<piece.cellData.size(); c++)
                {
                    printDoubleArray(f, piece.cellData[c].name, piece.cellData[c].data);
                }
                f << "\t\t\t</CellData>\n";
            }
            if (piece.pntData.size() != 0)
            {
                f << "\t\t\t<PointData Scalars=\"scalars\">\n";
                for (pieces_index_type a=0; a<piece.pntData.size(); a++)
                {
                    printDoubleArray(f, piece.pntData[a].name, piece.pntData[a].data);
                }
                f << "\t\t\t</PointData>\n";
            }
            f << "\t\t\t<Points>\n";
            printDoubleArray(f, "Position", piece.pnts);
            f << "\t\t\t</Points>\n";
            f << "\t\t\t<Cells>\n";
            printIntArray(f, "connectivity", piece.connectivity);
            Eigen::MatrixXi offset(piece.connectivity.rows(),1);
            Eigen::MatrixXi type(piece.connectivity.rows(),1);
            
            
            int verts = (int)piece.connectivity.cols();
            
            
            for (int i=0; i<piece.connectivity.rows(); i++)
            {
                offset(i) = (i+1)*verts;
                if (verts == 1)
                {
                    type(i) = 1;
                }
                else if (verts == 2)
                {
                    type(i) = 3;
                }
                else if (verts == 3)
                {
                    type(i) = 5;
                }
                else if (verts == 4)
                {
                    type(i) = 9;
                }
                else if (verts == 8)
                {
                    type(i) = 11;
                }
                else
                {
                    std::cout << "ERROR : Unknown cell type for VTU file format" << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
            printIntArray(f, "offsets", offset);
            printIntArray(f, "types", type);
            f << "\t\t\t</Cells>\n";
            f << "\t\t</Piece>\n";
        }
        f << "\t</UnstructuredGrid>\n";
        f << "</VTKFile>\n";
        f.close();
    }
    else
    {
        std::cout << "ERROR : Could not open VTU file" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void VTUfile::printDoubleArray(std::ofstream &f,std::string nname,Eigen::MatrixXd array)
{
    f << "\t\t\t\t<DataArray type=\"Float64\" Name=\"" << nname << "\" NumberOfComponents=\"" << array.cols() << "\" Format=\"ascii\">\n";
    for (int i=0; i<array.rows(); i++)
    {
        for (int j=0; j<array.cols(); j++)
        {
            f << array(i,j) << "\t";
        }
        f << "\n";
    }
    f << "\t\t\t\t</DataArray>\n";
}

void VTUfile::printIntArray(std::ofstream &f,std::string nname,Eigen::MatrixXi array)
{
    f << "\t\t\t\t<DataArray type=\"Int32\" Name=\"" << nname << "\" Format=\"ascii\">\n";
    for (int i=0; i<array.rows(); i++)
    {
        for (int j=0; j<array.cols(); j++)
        {
            f << array(i,j) << "\t";
        }
        f << "\n";
    }
    f << "\t\t\t\t</DataArray>\n";
}
