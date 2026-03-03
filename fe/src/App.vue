<template>
  <div id="app">
    <div class="Gallery">
      <ViewInstance
        v-for="drawData in drawDataList"
        :key="drawData.name"
        :drawData="drawData">
      </ViewInstance>
    </div>
  </div>
</template>

<script>
import *  as d3 from 'd3'
import axios from "axios";

import ViewInstance from '@/components/ViewInstance.vue'
export default {
  name: 'App',
  components: {
    ViewInstance
  },
  data(){
    return {
      drawDataList:[]
    }
  },
  methods:{
  },
  mounted(){
    // d3.json('static/Pancreatic_byGenes_20neighbors_2top_velo.json').then(res=>{
    //     let data = res;
    //     this.drawDataList.push(data)
    //   })
    
    axios({
      method:"post",
      url:"/api/init_plot_GridVis1",
      data:{
        'project_name':['celldancer_Pancreas','Deepvelo_Pancreas','scvelo_deterministic_Pancreas','scvelo_dynamic_Pancreas','unitvelo_Pancreas']
      }
    }).then(res=>{
      let data = res.data;
      console.log('data:',data)
      this.drawDataList = [data]
    }).catch(err=>{
      console.log('err:',err)
    })

  }
}
</script>

<style>
#app {
  font-family: 'Avenir', Helvetica, Arial, sans-serif;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  text-align: center;
  color: #2c3e50;
  width: 100%;
  height: 100%;
}

.gallery {
  display: flex;
  justify-content:space-around;
  flex-wrap: wrap;
}

</style>
