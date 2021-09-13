Shader "Unlit/triangleShader"
{
    SubShader
    {
        Tags { "RenderType" = "Opaque" }
        Cull Off

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag

            #include "UnityCG.cginc"

            struct appdata
            {
                uint vertexId : SV_VertexID;
                uint instanceId : SV_InstanceID;
            };

            struct v2f
            {
                float4 vertex : SV_POSITION;
                float3 normal : TEXCOORD0;
            };
            StructuredBuffer<float> position_buffer;
            StructuredBuffer<float> normal_buffer;
            StructuredBuffer<int> index_buffer;
            float3 getVertexPos(appdata v)
            {
                int idx = index_buffer[v.instanceId * 3 + v.vertexId] * 3;
                float3 ret;
                ret.x = position_buffer[idx + 0];
                ret.y = position_buffer[idx + 1];
                ret.z = position_buffer[idx + 2];
                return ret;
            }
            float3 getVertexNormal(appdata v)
            {
                float3 nor;
                int idx = index_buffer[v.instanceId * 3 + v.vertexId] * 3;
                nor.x = normal_buffer[idx + 0];
                nor.y = normal_buffer[idx + 1];
                nor.z = normal_buffer[idx + 2];
                return nor;
            }
            v2f vert(appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(getVertexPos(v));
                o.normal = UnityObjectToClipPos(getVertexNormal(v));
                //o.normal = getVertexNormal(v);
                return o;
            }

            fixed4 frag(v2f i) : SV_Target
            {
                float3 normal = i.normal;
                return fixed4(i.normal.z,0,0,1);
            }
            ENDCG
        }
    }
}
