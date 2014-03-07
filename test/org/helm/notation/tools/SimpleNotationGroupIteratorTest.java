package org.helm.notation.tools;

import static org.junit.Assert.*;

import java.util.ArrayList;

import org.junit.Test;

public class SimpleNotationGroupIteratorTest {

	@Test
	public void testSimplePeptide(){ 
		SimpleNotationGroupIterator iterator=new SimpleNotationGroupIterator("A.K.A");
		ArrayList<String> list=new ArrayList<String>();
		while (iterator.hasNextGroup()){
			String group=iterator.nextGroup();
			list.add(group);
			System.out.println(group);
			
		}
		assertEquals("there should be three groups!",3, list.size());
		assertEquals("A",list.get(0));
		assertEquals("K",list.get(1));
		assertEquals("A",list.get(2));		
		
		
     
		
		
	}
	
	@Test
	public void testSimpleRNA(){
		SimpleNotationGroupIterator iterator=new SimpleNotationGroupIterator("P.R(A)[sP].RP.R(G)P.[LR]([5meC])");
		
		ArrayList<String> list=new ArrayList<String>();
		while (iterator.hasNextGroup()){
			String group=iterator.nextGroup();
			list.add(group);
			System.out.println(group);
			
		}
		assertEquals("there should be five groups!",5, list.size());
		assertEquals("P",list.get(0));
		assertEquals("R(A)[sP]",list.get(1));
		assertEquals("RP",list.get(2));		
		assertEquals("R(G)P",list.get(3));
		assertEquals("[LR]([5meC])",list.get(4));
		
		

	}

	@Test
	public void testSimpleChem(){
		SimpleNotationGroupIterator iterator=new SimpleNotationGroupIterator("SMPEG2");
		String group=null;
		int i=0;
		while (iterator.hasNextGroup()){
			group=iterator.nextGroup();
			i++;
			
		}
		assertEquals("SMPEG2", group);
		assertEquals("there should be one group!",1, i);
	}
	
	@Test
	public void testInlineSmilesPeptide(){
		ArrayList<String> list=new ArrayList<String>();
		SimpleNotationGroupIterator iterator=new SimpleNotationGroupIterator("G.G.K.A.A.[[SeH]C[C@H](N[*])C([*])=O |$;;;;_R1;;_R2;$|].[meC]");
		while (iterator.hasNextGroup()){
			String group=iterator.nextGroup();
			list.add(group);
		}
		assertEquals("there should be 7 groups!",7, list.size());
		assertEquals("G",list.get(0));
		assertEquals("G",list.get(1));
		assertEquals("K",list.get(2));		
		assertEquals("A",list.get(3));
		assertEquals("A",list.get(4));
		assertEquals("[[SeH]C[C@H](N[*])C([*])=O |$;;;;_R1;;_R2;$|]",list.get(5));
		assertEquals("[meC]",list.get(6));
		
	}
	
	@Test
	public void testInlineSmilesChem(){

		SimpleNotationGroupIterator iterator=new SimpleNotationGroupIterator("C[C@H](N[*])C([*])=O |$;;;;_R1;;_R2;$|");
		String group=null;
		int i=0;
		while (iterator.hasNextGroup()){
			group=iterator.nextGroup();
			i++;
			
		}
		assertEquals("C[C@H](N[*])C([*])=O |$;;;;_R1;;_R2;$|", group);
		assertEquals("there should be one group!",1, i);
	}
		
	
	@Test
	public void testInlineSmilesRNA(){
		SimpleNotationGroupIterator iterator=new SimpleNotationGroupIterator("P.[[*]OC[C@@H]1CN([*])C[C@H]([*])O1 |$_R1;;;;;;_R2;;;_R3;$|](A)[sP].RP.R(G)P.[LR]([5meC])");
		
		ArrayList<String> list=new ArrayList<String>();
		while (iterator.hasNextGroup()){
			String group=iterator.nextGroup();
			list.add(group);
			System.out.println(group);
			
		}
		assertEquals("there should be five groups!",5, list.size());
		assertEquals("P",list.get(0));
		assertEquals("[[*]OC[C@@H]1CN([*])C[C@H]([*])O1 |$_R1;;;;;;_R2;;;_R3;$|](A)[sP]",list.get(1));
		assertEquals("RP",list.get(2));		
		assertEquals("R(G)P",list.get(3));
		assertEquals("[LR]([5meC])",list.get(4));
		
		

	}
	
	
	
	
	

}
