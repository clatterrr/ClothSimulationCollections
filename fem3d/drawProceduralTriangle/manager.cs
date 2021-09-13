using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;

public class manager : MonoBehaviour
{
    public Material material;
    float[] position_data = new float[18];

    const int Nx = 2;
    const int Ny = Nx;
    const int node_num = Nx * Nx;
    const int element_num = (Nx - 1) * (Nx - 1) * 2;

    float[] node_pos_data = new float[node_num * 3];
    int[] element_data = new int[element_num * 3];


    ComputeBuffer node_pos;
    ComputeBuffer element;
    private void Start()
    {
        node_pos = new ComputeBuffer(node_num * 3, sizeof(float));
        element = new ComputeBuffer(element_num * 3, sizeof(int));

        int cnt = 0;
        for(int j = 0;j  < Ny;j++)
        {
            for(int i = 0;i < Nx;i++)
            {
                node_pos_data[cnt + 0] = i;
                node_pos_data[cnt + 1] = j;
                node_pos_data[cnt + 2] = 0;
                cnt += 3;
            }
        }
        node_pos.SetData(node_pos_data);
        cnt = 0;
        for (int j = 0; j < (Ny - 1); j++)
        {
            for (int i = 0; i < (Nx - 1); i++)
            {
                int idx0 = j * Nx + i;
                int idx1 = idx0 + 1;
                int idx2 = idx0 + Nx;
                int idx3 = idx2 + 1;
                element_data[cnt + 0] = idx0;
                element_data[cnt + 1] = idx1;
                element_data[cnt + 2] = idx2;

                element_data[cnt + 3] = idx2;
                element_data[cnt + 4] = idx1;
                element_data[cnt + 5] = idx3;
                cnt += 6;
            }
        }
        element.SetData(element_data);
    }

    private void Update()
    {

        
    }
    void OnRenderObject()
    {
        //https://zenn.dev/fuqunaga/scraps/782ea5ca0b002f
        // 花了四个小时，终于找到一个简单的教程了
        material.SetBuffer("position_buffer", node_pos);
        material.SetBuffer("index_buffer", element);
        material.SetPass(0);
        Graphics.DrawProceduralNow(MeshTopology.Triangles, 3, element_num);
    }
}
